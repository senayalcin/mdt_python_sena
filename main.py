import argparse
import gzip
from collections import defaultdict
import sys

import args
from scipy.stats import binomtest



def parse_arguments():
    parser = argparse.ArgumentParser(
        description="This script analyzes dominant transcripts by comparing PCAWG and GTEx RNA-seq data.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-pcawg", required=True, help="PCAWG RNA-seq TPM file (e.g., pcawg.rnaseq.transcript.expr.tpm.tsv)")
    parser.add_argument("-gtex", required=True, help="GTEx RNA-seq TPM file")
    parser.add_argument("-ensgFile", required=True, help="ENSG-ENSP-ENST mapping file")
    parser.add_argument("-canonFile", required=True, help="Canonical isoform information")
    parser.add_argument("-isoformIntFile", required=True, help="Isoform interaction file")
    parser.add_argument("-sequenceFile", required=True, help="ENSP sequence file")
    parser.add_argument("-minStringScore", type=int, required=True, help="Minimum STRING score (e.g., 900)")
    parser.add_argument("-enstColumnNuo", type=int, required=True, help="ENST column index in the TPM file (e.g., 0)")
    parser.add_argument("-minEnrichment", type=float, required=True, help="Minimum enrichment factor (e.g., 2.0)")
    parser.add_argument("-maxQvalue", type=float, required=True, help="Maximum q-value (e.g., 0.01)")
    parser.add_argument("-minExpress", type=float, required=True, help="Expression threshold (e.g., 2.0)")

    parser.add_argument("-pvalueFile", help="Optional p-value file")
    parser.add_argument("-randomPvalueFile", help="Random p-value file")
    parser.add_argument("-pcawgEnrichmentFile", help="PCAWG enrichment output")
    parser.add_argument("-gtexEnrichmentFile", help="GTEx enrichment output")
    parser.add_argument("-redundantFile", help="Redundant ENST file")
    parser.add_argument("-verbose", action="store_true", help="Enable verbose output")

    return parser.parse_args()



def read_enst2ensg_file(enst2ensg_filepath):
    """
    Returns a dictionary containing ENST → ENSG mappings.
    The file is in gzip format and tab-delimited.
    """
    enst2ensg = {}

    with gzip.open(enst2ensg_filepath, "rt") as f:
        for line in f:
            if not line.startswith("ENS"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                enst_id = parts[2]
                ensg_id = parts[0]
                enst2ensg[enst_id] = ensg_id

    return enst2ensg



def read_expression_file(file_path, enst2ensg, redundant_dict, enst_col_index=0, verbose=False):
    """
    Reads a TPM file and returns:
    - gene-level expression
    - transcript-level expression
    - relative expressions
    """

    enst_express = defaultdict(lambda: defaultdict(dict))
    ensg_express = defaultdict(lambda: defaultdict(float))

    with gzip.open(file_path, "rt") as f:
        lines = f.readlines()

    header = lines[0].strip().split("\t")
    samples = header

    for line in lines[1:]:
        parts = line.strip().split("\t")
        if len(parts) <= enst_col_index:
            continue

        enst = parts[enst_col_index]

        # Redundant correction
        if enst in redundant_dict:
            if verbose:
                print(f"Replacing redundant cDNA {enst} with {redundant_dict[enst]}")
            enst = redundant_dict[enst]

        # Format check
        if "ENST0" not in enst:
            continue

        enst = enst.split(".")[0]  # ENST00000123456.1 → ENST00000123456

        if enst not in enst2ensg:
            if verbose:
                print(f"WARNING: {enst} not in enst2ensg mapping file. Skipping.")
            continue

        ensg = enst2ensg[enst]

        for i in range(enst_col_index + 1, len(parts)):
            expression = parts[i]
            expression = 0.0 if expression == "NA" else float(expression)

            sample = samples[i]

            enst_express[ensg][sample][enst] = enst_express[ensg][sample].get(enst, 0) + expression
            ensg_express[ensg][sample] += expression

    rel_enst_express = defaultdict(lambda: defaultdict(dict))
    rel_enst_express_string = defaultdict(lambda: defaultdict(str))
    all_rel_enst_express = []

    for ensg in enst_express:
        for sample in enst_express[ensg]:
            gene_expr = ensg_express[ensg][sample]
            for enst in enst_express[ensg][sample]:
                trans_expr = enst_express[ensg][sample][enst]
                rel_expr = trans_expr / gene_expr if gene_expr > 0 else 0.0

                rel_enst_express[ensg][enst][sample] = rel_expr
                all_rel_enst_express.append(rel_expr)

                if rel_enst_express_string[ensg][enst]:
                    rel_enst_express_string[ensg][enst] += f" {rel_expr:.3f}"
                else:
                    rel_enst_express_string[ensg][enst] = f"{rel_expr:.3f}"

    return ensg_express, enst_express, rel_enst_express, rel_enst_express_string, all_rel_enst_express



def compute_enrichment(expression_data, min_expression, min_enrichment, output_file_path=None):
    """
    Calculates dominant transcript enrichment ratio.

    Parameters:
    - expression_data: dict[ensg][sample][enst] = expression (TPM)
    - min_expression: float, minimum MDT expression (e.g., 2.0)
    - min_enrichment: float, minimum fold change required for MDT over the second transcript
    - output_file_path: str, if provided, writes results to a gzip-compressed .tsv.gz file

    Returns:
    - enrichment_dict: dict[ensg][dominant_enst][sample] = enrichment value
    - all_enrichment_values: list[float], all enrichment values for statistical testing
    """

    enrichment_dict = defaultdict(lambda: defaultdict(dict))
    all_enrichment_values = []

    output_file = gzip.open(output_file_path, "wt") if output_file_path else None

    for ensg in sorted(expression_data):
        for sample in sorted(expression_data[ensg]):

            ensts = expression_data[ensg][sample]
            sorted_ensts = sorted(ensts.items(), key=lambda x: x[1], reverse=True)

            if len(sorted_ensts) < 2:
                continue

            # dominant and second transcript
            dominant_enst, dominant_expr = sorted_ensts[0]
            other_enst, other_expr = sorted_ensts[1]

            # Is dominant expression above the threshold?
            if dominant_expr >= min_expression:
                enrichment = dominant_expr / other_expr if other_expr > 0 else dominant_expr

                if output_file:
                    output_file.write(f"{sample}\t{ensg}\t{dominant_enst}\t{other_enst}\t{dominant_expr:.3f}\t{other_expr:.3f}\t{enrichment:.3f}\n")

                if enrichment >= min_enrichment:
                    enrichment_dict[ensg][dominant_enst][sample] = enrichment
                    all_enrichment_values.append(enrichment)

    if output_file:
        output_file.close()

    return enrichment_dict, all_enrichment_values


def read_canon_file(canon_file_path):
    """
    Reads canonical protein IDs. Performs three types of mappings:
    - ENST → Canonical ENSP
    - ENSG → Canonical ENSP
    - ENSP → ENSG
    """

    canon = {}

    with gzip.open(canon_file_path, "rt") as f:
        for line in f:
            if not line.startswith("ENS"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue

            canon_ensp, ensg, ensp, enst, gene_name, source = parts[:6]

            canon[enst] = canon_ensp
            canon[ensg] = canon_ensp
            canon[ensp] = ensg

    return canon


def read_isoform_interaction_file(isoform_file_path, min_string_score):
    """
    For each transcript, returns:
    - Missing interactions (stringScore ≥ min)
    - Total number of interactions (missing + existing)

    Returns:
    - miss_int_for_isoform: dict[ENST][interaction_str] = 1
    - interaction_count: dict[ENST] = total count
    """

    miss_int_for_isoform = defaultdict(dict)
    interaction_count = defaultdict(int)

    with gzip.open(isoform_file_path, "rt") as f:
        for line in f:
            if not line.startswith("EN"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue

            enst1 = parts[2]
            miss_int_field = parts[9]
            exist_int_field = parts[8]

            interaction_count[enst1] = 0

            # Missing interactions
            miss_interactions = miss_int_field.split(",")
            for int_str in miss_interactions:
                if ":" in int_str:
                    fields = int_str.split(":")
                    if len(fields) >= 8:
                        string_score = float(fields[6])
                        if string_score >= min_string_score:
                            miss_int_for_isoform[enst1][int_str] = 1
                            interaction_count[enst1] += 1

            # Existing interactions
            exist_interactions = exist_int_field.split(",")
            for int_str in exist_interactions:
                if ":" in int_str:
                    fields = int_str.split(":")
                    if len(fields) >= 8:
                        interaction_count[enst1] += 1

    return miss_int_for_isoform, interaction_count

def read_redundant_file(redundant_file_path, verbose=False):
    """
    Maps redundant transcripts.
    enst_redundant → enst_selected
    """
    redundant_map = {}

    with open(redundant_file_path, "r") as f:
        for line in f:
            if not line.startswith("EN"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                ensg, enst_sel, enst_red = parts[:3]
                redundant_map[enst_red] = enst_sel

    if verbose:
        print(f"Found {len(redundant_map)} redundant ENSTs", file=sys.stderr)

    return redundant_map


def read_sequence_file(sequence_file_path, verbose=False):
    """
    Returns a mapping of ENST → protein/amino acid sequence.
    """
    sequences = {}

    with gzip.open(sequence_file_path, "rt") as f:
        for line in f:
            if not line.startswith("EN"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                ensp, ensg, enst, sequence = parts[:4]
                sequences[enst] = sequence

    if verbose:
        print(f"Found {len(sequences)} sequences", file=sys.stderr)

    return sequences

def median(values):
    """
    Returns the median of a list of numbers.
    """
    sorted_vals = sorted(values)
    n = len(sorted_vals)

    if n == 0:
        return 0  # or None?

    if n % 2 == 1:
        return sorted_vals[n // 2]
    else:
        return (sorted_vals[(n // 2) - 1] + sorted_vals[n // 2]) / 2



def main():
    args = parse_arguments()

    # Check for required parameters (argparse already does this, but we include it for safety)
    required_params = [
        args.pcawg, args.gtex, args.ensgFile, args.canonFile,
        args.isoformIntFile, args.sequenceFile, args.enstColumnNuo,
        args.minStringScore, args.minEnrichment, args.maxQvalue, args.minExpress
    ]
    if not all(required_params):
        print("Missing required parameters. Please use -h for help.", file=sys.stderr)
        sys.exit(1)

    # Redundant and mapping files
    enst2ensg = read_enst2ensg_file(args.ensgFile)
    redundant_map = read_redundant_file(args.redundantFile, verbose=args.verbose)
    sequences = read_sequence_file(args.sequenceFile, verbose=args.verbose)
    canonical = read_canon_file(args.canonFile)

    # Expression files
    if args.verbose:
        print(f"Reading PCAWG expression file: {args.pcawg}", file=sys.stderr)
    ensg_pcawg, enst_pcawg, rel_pcawg, relstr_pcawg, allrel_pcawg = read_expression_file(
        args.pcawg, enst2ensg=enst2ensg, redundant_dict=redundant_map,
        enst_col_index=args.enstColumnNuo, verbose=args.verbose
    )

    if args.verbose:
        print(f"Reading GTEx expression file: {args.gtex}", file=sys.stderr)
    ensg_gtex, enst_gtex, rel_gtex, relstr_gtex, allrel_gtex = read_expression_file(
        args.gtex, enst2ensg=enst2ensg, redundant_dict=redundant_map,
        enst_col_index=args.enstColumnNuo, verbose=args.verbose
    )

    # Interaction file
    miss_iso, int_counts = read_isoform_interaction_file(args.isoformIntFile, args.minStringScore)

    # Enrichment computation
    enrichment_pcawg, all_enrich_pcawg = compute_enrichment(
        enst_pcawg, args.minExpress, args.minEnrichment, output_file_path=args.pcawgEnrichmentFile
    )

    enrichment_gtex, all_enrich_gtex = compute_enrichment(
        enst_gtex, args.minExpress * 0.1, args.minEnrichment, output_file_path=args.gtexEnrichmentFile
    )

    header = (
    "#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\t"
    "GTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\tNumberOfGTExMDT\tcMDTenrichment\t"
    "GTExMDTmedianRelExpression\tcMDTrelExpression\t"
    "cMDTuniqMissedInt\tTotalNumberOfSTRINGint\t"
    "cMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt"
    )
    print(header)


    pvalue_file = None
    random_pvalue_file = None

    if args.pvalueFile:
        pvalue_file = open(args.pvalueFile, "w")
        pvalue_file.write(
            "#ENSG\tcMDT\tRNAseqAliquotID\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
            "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\tRelMDTexpressionInGTEx\n"
        )

    if args.randomPvalueFile:
        random_pvalue_file = open(args.randomPvalueFile, "w")
        random_pvalue_file.write(
            "#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
            "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n"
        )




