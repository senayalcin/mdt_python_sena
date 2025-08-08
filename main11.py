#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of mostDominantTranscriptSwitch.pl
Author: Python conversion
Date: 2025

Determines cancer-specific MDT based on 2-fold expression difference to
minor dominant transcripts and subsequently compares the relative
expression difference to expression differences in GTEx using Sign-test.
"""

import argparse
import gzip
import sys
import numpy as np
import random
from scipy import stats
from scipy.stats import binomtest  # For newer scipy versions
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings('ignore')

# Set random seed for reproducibility (Perl uses srand(23))
random.seed(23)
np.random.seed(23)

# Global variables
enst_column_no = None


def print_help():
    """Print help message"""
    help_text = """
Usage: python mostDominantTranscriptSwitch.py [options]

Options:
    -h, --help                      Show this help message and exit
    --pcawg FILE                    PCAWG expression file (required)
    --gtex FILE                     GTEx expression file (required)
    --ensg FILE                     ENSG to ENST mapping file (required)  
    --canon FILE                    Canonical isoforms file (required)
    --iso FILE                      Isoform interaction file (required)
    --seq FILE                      Sequence file (required)
    --redundant FILE                Redundant ENST file (optional)
    --minString INT                 Minimum STRING score (required)
    --enst INT                      ENST column number (required)
    --minEnr FLOAT                  Minimum enrichment fold (required)
    --maxQ FLOAT                    Maximum Q-value (required)
    --minExp FLOAT                  Minimum expression value (required)
    --pvalueFile FILE               P-value output file (optional)
    --randomPvalueFile FILE         Random p-value output file (optional)
    --pcawgEnrichmentFile FILE      PCAWG enrichment output file (optional)
    --gtexEnrichmentFile FILE       GTEx enrichment output file (optional)
    --verbose                       Verbose output
    """
    print(help_text)


def read_enst2ensg_file(filename):
    """Read ENST to ENSG mapping file"""
    enst2ensg = {}

    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'):
                continue
            parts = line.strip().split('\t')
            if len(parts) > 2:
                enst2ensg[parts[2]] = parts[0]

    return enst2ensg


def read_redundant_file(filename, verbose=False):
    """Read redundant ENST file"""
    redundant = {}

    if not filename:
        return redundant

    try:
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('EN'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    ensg, enst_sel, enst_red = parts[:3]
                    redundant[enst_red] = enst_sel
    except:
        pass

    if verbose and redundant:
        print(f"Found {len(redundant)} redundant ENSTs", file=sys.stderr)

    return redundant


def read_sequence_file(filename, verbose=False):
    """Read sequence file"""
    sequences = {}

    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                ensp, ensg, enst, seq = parts[:4]
                sequences[enst] = seq

    if verbose:
        print(f"Found {len(sequences)} sequences", file=sys.stderr)

    return sequences


def read_expression_file(filename, sequences, enst2ensg, redundant, verbose=False):
    """Read expression file and calculate relative expressions"""
    enst_express = {}
    ensg_express = {}

    with gzip.open(filename, 'rt') as f:
        # Read header
        header = f.readline().strip()
        samples = header.split('\t')

        for line in f:
            parts = line.strip().split('\t')

            # Skip empty lines or lines with insufficient columns
            if len(parts) <= enst_column_no:
                continue

            enst = parts[enst_column_no]

            # Handle redundant transcripts
            if enst in redundant:
                if verbose:
                    print(f"Replacing redundant cDNA {enst} with {redundant[enst]}")
                enst = redundant[enst]

            if 'ENST0' not in enst:
                continue

            # Remove version number
            enst = enst.split('.')[0]

            if enst not in enst2ensg:
                if verbose:
                    print(f"WARNING: {enst} does not exist in ENST2ENSG file. Skipping transcript.", file=sys.stderr)
                continue

            ensg = enst2ensg[enst]

            # Process expression values for each sample
            # Make sure we don't exceed the bounds of both lists
            for i in range(enst_column_no + 1, min(len(parts), len(samples))):
                expression = parts[i]
                if expression == "NA":
                    expression = 0.0
                else:
                    try:
                        expression = float(expression)
                    except ValueError:
                        expression = 0.0

                sample = samples[i]

                # Initialize dictionaries if needed
                if ensg not in enst_express:
                    enst_express[ensg] = {}
                if sample not in enst_express[ensg]:
                    enst_express[ensg][sample] = {}
                if ensg not in ensg_express:
                    ensg_express[ensg] = {}

                # Add expression values
                if enst in enst_express[ensg][sample]:
                    enst_express[ensg][sample][enst] += expression
                else:
                    enst_express[ensg][sample][enst] = expression

                if sample in ensg_express[ensg]:
                    ensg_express[ensg][sample] += expression
                else:
                    ensg_express[ensg][sample] = expression

    # Calculate relative transcript expression
    rel_enst_express = {}
    rel_enst_express_string = {}
    all_rel_enst_express = []

    for ensg in sorted(enst_express.keys()):
        for sample in sorted(enst_express[ensg].keys()):
            gene_express = ensg_express[ensg][sample]

            for enst in sorted(enst_express[ensg][sample].keys()):
                trans_express = enst_express[ensg][sample][enst]

                rel_express = 0.0
                if gene_express > 0:
                    rel_express = trans_express / gene_express

                if ensg not in rel_enst_express:
                    rel_enst_express[ensg] = {}
                if enst not in rel_enst_express[ensg]:
                    rel_enst_express[ensg][enst] = {}
                rel_enst_express[ensg][enst][sample] = rel_express
                all_rel_enst_express.append(rel_express)

                if ensg not in rel_enst_express_string:
                    rel_enst_express_string[ensg] = {}
                if enst not in rel_enst_express_string[ensg]:
                    rel_enst_express_string[ensg][enst] = ""

                # Build space-separated string of relative expressions
                if rel_enst_express_string[ensg][enst]:
                    rel_enst_express_string[ensg][enst] += f" {rel_express:.3f}"
                else:
                    rel_enst_express_string[ensg][enst] = f"{rel_express:.3f}"

    return ensg_express, enst_express, rel_enst_express, rel_enst_express_string, all_rel_enst_express


def enrichment(express, minimum_express, min_enrichment, output_file=None):
    """Calculate enrichment for dominant transcripts"""
    enrichment_dict = {}
    all_enrichment = []

    out_f = None
    if output_file:
        out_f = gzip.open(output_file, 'wt')

    for ensg in sorted(express.keys()):
        for sample in sorted(express[ensg].keys()):
            most_dominant_enst = "-"
            most_dominant_enst_express = 0

            # Sort isoforms by highest expression
            sorted_enst = sorted(express[ensg][sample].keys(),
                                 key=lambda x: express[ensg][sample][x],
                                 reverse=True)

            n = 0
            for enst in sorted_enst:
                enst_express = express[ensg][sample][enst]

                # First transcript is potentially the MDT
                n += 1
                if n == 1:
                    most_dominant_enst = enst
                    most_dominant_enst_express = express[ensg][sample][enst]
                    continue

                enrichment_val = most_dominant_enst_express
                # Ignore transcripts with insignificant expression
                if most_dominant_enst_express >= minimum_express:
                    if enst_express > 0:
                        enrichment_val = most_dominant_enst_express / enst_express

                    if out_f:
                        out_f.write(f"{sample}\t{ensg}\t{most_dominant_enst}\t{enst}\t"
                                    f"{most_dominant_enst_express}\t{enst_express}\t{enrichment_val}\n")

                    if enrichment_val >= min_enrichment:
                        if ensg not in enrichment_dict:
                            enrichment_dict[ensg] = {}
                        if most_dominant_enst not in enrichment_dict[ensg]:
                            enrichment_dict[ensg][most_dominant_enst] = {}
                        enrichment_dict[ensg][most_dominant_enst][sample] = enrichment_val
                        all_enrichment.append(enrichment_val)

                # Only compare first two transcripts
                if n == 2:
                    break

    if out_f:
        out_f.close()

    return enrichment_dict, all_enrichment


def read_canon_file(filename):
    """Read canonical isoforms file"""
    canon = {}

    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                canon_ensp, ensg, ensp, enst, gene_name, source = parts[:6]
                canon[enst] = canon_ensp
                canon[ensg] = canon_ensp
                canon[ensp] = ensg

    return canon


def read_isoform_int_file(filename, min_string_score):
    """Read isoform interaction file"""
    miss_int4isoform = {}
    int_n = {}

    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                ensg1, ensp1, enst1, string_ensp1, string_int_n, exist_int_n, miss_int_n, rel_miss_int_n, exist_int, miss_int = parts[
                                                                                                                                :10]
                int_n[enst1] = 0

                if enst1 not in miss_int4isoform:
                    miss_int4isoform[enst1] = {}

                # Process missed interactions
                if miss_int and miss_int != "":
                    for interaction in miss_int.split(','):
                        if ':' in interaction:
                            int_parts = interaction.split(':')
                            if len(int_parts) >= 8:
                                try:
                                    string_score = float(int_parts[6])
                                    if string_score >= min_string_score:
                                        miss_int4isoform[enst1][interaction] = 1
                                        int_n[enst1] += 1
                                except:
                                    pass

                # Process existing interactions
                if exist_int and exist_int != "":
                    for interaction in exist_int.split(','):
                        if ':' in interaction:
                            int_parts = interaction.split(':')
                            if len(int_parts) >= 8:
                                try:
                                    string_score = float(int_parts[6])
                                    if string_score >= min_string_score:
                                        int_n[enst1] += 1
                                except:
                                    pass

    return miss_int4isoform, int_n


def median(values):
    """Calculate median of a list - matching Perl implementation"""
    if not values:
        return 0

    sorted_values = sorted(values)
    n = len(sorted_values)

    # Perl uses different logic: if odd number, take middle; if even, average of two middle
    if n % 2 == 1:  # Odd number of elements
        return sorted_values[n // 2]
    else:  # Even number of elements
        return (sorted_values[n // 2 - 1] + sorted_values[n // 2]) / 2


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Dominant Transcript Switch Analysis')
    parser.add_argument('--pcawg', required=True, help='PCAWG expression file')
    parser.add_argument('--gtex', required=True, help='GTEx expression file')
    parser.add_argument('--ensg', required=True, help='ENSG to ENST mapping file')
    parser.add_argument('--canon', required=True, help='Canonical isoforms file')
    parser.add_argument('--iso', required=True, help='Isoform interaction file')
    parser.add_argument('--seq', required=True, help='Sequence file')
    parser.add_argument('--redundant', help='Redundant ENST file')
    parser.add_argument('--minString', type=int, required=True, help='Minimum STRING score')
    parser.add_argument('--enst', type=int, required=True, help='ENST column number')
    parser.add_argument('--minEnr', type=float, required=True, help='Minimum enrichment fold')
    parser.add_argument('--maxQ', type=float, required=True, help='Maximum Q-value')
    parser.add_argument('--minExp', type=float, required=True, help='Minimum expression value')
    parser.add_argument('--pvalueFile', help='P-value output file')
    parser.add_argument('--randomPvalueFile', help='Random p-value output file')
    parser.add_argument('--pcawgEnrichmentFile', help='PCAWG enrichment output file')
    parser.add_argument('--gtexEnrichmentFile', help='GTEx enrichment output file')
    parser.add_argument('--verbose', action='store_true', help='Verbose output')

    args = parser.parse_args()

    global enst_column_no
    enst_column_no = args.enst

    # Read files
    if args.verbose:
        print(f"Reading {args.seq} ... ", end='', file=sys.stderr)
    sequences = read_sequence_file(args.seq, args.verbose)
    if args.verbose:
        print("done", file=sys.stderr)

    # Read ENST2ENSG mapping
    enst2ensg = read_enst2ensg_file(args.ensg)

    # Read redundant file
    redundant = read_redundant_file(args.redundant, args.verbose) if args.redundant else {}

    if args.verbose:
        print(f"Reading in {args.pcawg} ... ", end='', file=sys.stderr)
    ensg_express_pcawg, enst_express_pcawg, rel_enst_express_pcawg, rel_enst_express_string_pcawg, all_rel_enst_express_pcawg = \
        read_expression_file(args.pcawg, sequences, enst2ensg, redundant, args.verbose)
    if args.verbose:
        print("done", file=sys.stderr)

    if args.verbose:
        print(f"Reading in {args.gtex} ... ", end='', file=sys.stderr)
    ensg_express_gtex, enst_express_gtex, rel_enst_express_gtex, rel_enst_express_string_gtex, all_rel_enst_express_gtex = \
        read_expression_file(args.gtex, sequences, enst2ensg, redundant, args.verbose)
    if args.verbose:
        print("done", file=sys.stderr)

    if args.verbose:
        print(f"Reading {args.canon} ... ", end='', file=sys.stderr)
    canonical = read_canon_file(args.canon)
    if args.verbose:
        print("done", file=sys.stderr)

    if args.verbose:
        print(f"Reading {args.iso} ... ", end='', file=sys.stderr)
    miss_iso_int, int_n = read_isoform_int_file(args.iso, args.minString)
    if args.verbose:
        print("done", file=sys.stderr)

    if args.verbose:
        print("Calculating Enrichment for PCAWG ... ", end='', file=sys.stderr)
    enrichment_pcawg, all_enrichment_pcawg = enrichment(enst_express_pcawg, args.minExp, args.minEnr,
                                                        args.pcawgEnrichmentFile)
    if args.verbose:
        print("done", file=sys.stderr)

    if args.verbose:
        print("Calculating Enrichment for GTEx ... ", end='', file=sys.stderr)
    enrichment_gtex, all_enrichment_gtex = enrichment(enst_express_gtex, args.minExp * 0.1, args.minEnr,
                                                      args.gtexEnrichmentFile)
    if args.verbose:
        print("done", file=sys.stderr)

    # Print header
    print("#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\t"
          "GTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\tNumberOfGTExMDT\tcMDTenrichment\t"
          "GTExMDTmedianRelExpression\tcMDTrelExpression\t"
          "cMDTuniqMissedInt\tTotalNumberOfSTRINGint\t"
          "cMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt")

    # Open output files if specified
    pvalue_file = None
    random_pvalue_file = None

    if args.pvalueFile:
        pvalue_file = open(args.pvalueFile, 'w')
        pvalue_file.write("#ENSG\tcMDT\ttRNAseqAliquotID\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
                          "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\tRelMDTexpressionInGTEx\n")

    if args.randomPvalueFile:
        random_pvalue_file = open(args.randomPvalueFile, 'w')
        random_pvalue_file.write(
            "#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
            "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n")

    # Process results
    output_part1 = []
    output_part2 = []
    pvalues_list = []

    for ensg in sorted(enrichment_pcawg.keys()):
        if ensg not in rel_enst_express_string_gtex:
            if args.verbose:
                print(f"WARNING: {ensg} does not exists in {args.gtex}. Skipping {ensg}.", file=sys.stderr)
            continue

        for enst in sorted(enrichment_pcawg[ensg].keys()):
            if enst not in rel_enst_express_string_gtex[ensg]:
                if args.verbose:
                    print(f"WARNING: {ensg} - {enst} does not exists in {args.gtex}. Skipping it.", file=sys.stderr)
                continue

            for sample in sorted(enrichment_pcawg[ensg][enst].keys()):
                enrichment_pcawg_val = enrichment_pcawg[ensg][enst][sample]
                enst_rel_express_pcawg = rel_enst_express_pcawg[ensg][enst][sample]
                enst_rel_expresses_gtex = rel_enst_express_string_gtex[ensg][enst]

                # Only do analysis for MDT that are specific to PCAWG
                if enst not in enrichment_gtex.get(ensg, {}) and ensg in enrichment_gtex and len(
                        enrichment_gtex[ensg]) > 0:
                    enst_rel_expresses_gtex_list = [float(x) for x in enst_rel_expresses_gtex.split()]

                    rel_expression_gtex_median = 0
                    if enst_rel_expresses_gtex_list:
                        rel_expression_gtex_median = median(enst_rel_expresses_gtex_list)

                    pos = sum(1 for x in enst_rel_expresses_gtex_list if enst_rel_express_pcawg > x)
                    neg = sum(1 for x in enst_rel_expresses_gtex_list if enst_rel_express_pcawg < x)

                    # Binomial test
                    if pos + neg > 0:
                        result = binomtest(pos, pos + neg, p=0.5, alternative='two-sided')
                        pvalue = result.pvalue
                    else:
                        pvalue = 1.0

                    if pvalue_file:
                        pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t"
                                          f"{enst_rel_express_pcawg:.3f}\t{pos}\t{neg}\t{pvalue:.3e}\t"
                                          f"{enst_rel_expresses_gtex}\n")

                    # Random p-value calculation (100 iterations like in Perl)
                    if random_pvalue_file:
                        for j in range(100):
                            random_enrichment_pcawg = random.choice(all_enrichment_pcawg)

                            r_pos = sum(1 for x in enst_rel_expresses_gtex_list if random_enrichment_pcawg > x)
                            r_neg = sum(1 for x in enst_rel_expresses_gtex_list if random_enrichment_pcawg < x)

                            if r_pos + r_neg > 0:
                                r_result = binomtest(r_pos, r_pos + r_neg, p=0.5, alternative='two-sided')
                                random_pvalue = r_result.pvalue
                            else:
                                random_pvalue = 1.0

                            random_pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t"
                                                     f"{random_enrichment_pcawg:.3f}\t{r_pos}\t{r_neg}\t"
                                                     f"{random_pvalue:.3e}\n")

                    # Only continue if p-value is significant (less than maxQvalue for initial filtering)
                    if pvalue < args.maxQ:
                        pvalues_list.append(pvalue)

                        # Check GTEx MDTs
                        gtex_mdts = {}
                        mdi_enrichments_gtex = []

                        if ensg in enrichment_gtex:
                            for n_enst in sorted(enrichment_gtex[ensg].keys()):
                                for n_sample in sorted(enrichment_gtex[ensg][n_enst].keys()):
                                    n_enrichment_gtex = enrichment_gtex[ensg][n_enst][n_sample]
                                    # Check if sequences are different
                                    if (n_enrichment_gtex >= args.minEnr and
                                            n_enst in sequences and enst in sequences and
                                            sequences[n_enst] != sequences[enst]):
                                        if n_enst not in gtex_mdts:
                                            gtex_mdts[n_enst] = 0
                                        gtex_mdts[n_enst] += 1
                                        mdi_enrichments_gtex.append(n_enrichment_gtex)

                        # Skip if GTEx doesn't have MDT for 50% of its cases
                        gtex_mdts_str = ""
                        skip = True
                        for g in sorted(gtex_mdts.keys(), key=lambda x: gtex_mdts[x], reverse=True):
                            gtex_mdts_str += f"{g}:{gtex_mdts[g]},"
                            if gtex_mdts[g] >= len(enst_rel_expresses_gtex_list) * 0.5:
                                skip = False

                        if skip:
                            continue

                        gtex_mdts_str = gtex_mdts_str.rstrip(',') if gtex_mdts_str else "-"

                        median_mdt_enrichment_gtex = "-"
                        if mdi_enrichments_gtex:
                            median_mdt_enrichment_gtex = f"{median(mdi_enrichments_gtex):.3f}"

                        # Check disrupted interactions
                        string_int_n = int_n.get(enst, -1)
                        missed_int = ""
                        missed_cancer_int_n = -1
                        missed_common_int_n = -1

                        if enst in miss_iso_int:
                            missed_cancer_int_n = 0
                            missed_common_int_n = 0

                            for c_int in sorted(miss_iso_int[enst].keys()):
                                int_parts = c_int.split(':')
                                if len(int_parts) >= 8:
                                    ensp2 = int_parts[5]

                                    # Check if interaction partner is expressed
                                    if ensp2 in canonical and canonical[ensp2] in ensg_express_pcawg:
                                        if sample in ensg_express_pcawg[canonical[ensp2]]:
                                            found = False

                                            # Check if interaction is disrupted in GTEx MDTs
                                            if ensg in enrichment_gtex:
                                                for n_enst in enrichment_gtex[ensg]:
                                                    for n_sample in enrichment_gtex[ensg][n_enst]:
                                                        enrichment_gtex_val = enrichment_gtex[ensg][n_enst][n_sample]
                                                        if enrichment_gtex_val >= args.minEnr:
                                                            if n_enst in miss_iso_int and c_int in miss_iso_int[n_enst]:
                                                                found = True
                                                                break
                                                    if found:
                                                        break

                                            if found:
                                                missed_common_int_n += 1
                                            else:
                                                missed_int += f"{c_int},"
                                                missed_cancer_int_n += 1

                        missed_int = missed_int.rstrip(',') if missed_int else "-"

                        # Store output
                        canon_ensp = canonical.get(enst, "-")

                        output_part1.append(f"{ensg}\t{canon_ensp}\t{enst}\t{sample}\t"
                                            f"{gtex_mdts_str}\t{len(enst_rel_expresses_gtex_list)}\t{median_mdt_enrichment_gtex}")

                        output_part2.append(f"{len(mdi_enrichments_gtex)}\t{enrichment_pcawg_val:.3f}\t"
                                            f"{rel_expression_gtex_median:.3f}\t{enst_rel_express_pcawg:.3f}\t"
                                            f"{missed_int}\t{string_int_n}\t"
                                            f"{missed_cancer_int_n}\t{missed_common_int_n}")

    # FDR correction
    if pvalues_list:
        # Use FDR correction
        reject, qvalues, alpha_sidak, alpha_bonf = multipletests(pvalues_list, method='fdr_bh')
    else:
        qvalues = []

    # Output results
    for i in range(len(output_part1)):
        if i < len(qvalues) and qvalues[i] < args.maxQ:
            print(f"{output_part1[i]}\t{pvalues_list[i]:.3e}\t{qvalues[i]:.3e}\t{output_part2[i]}")

    # Close files
    if pvalue_file:
        pvalue_file.close()
    if random_pvalue_file:
        random_pvalue_file.close()


if __name__ == "__main__":
    main()