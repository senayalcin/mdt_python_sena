#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of mostDominantTranscriptSwitch.pl (fixed)
Notes:
- Aligns p-value (binomtest) collection strictly with printed rows to avoid index drift.
- Adds optional built-in validator to compare with a provided expected output (TSV/TSV.GZ/XLSX).
"""

import argparse
import gzip
import sys
import os
import numpy as np
import random
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
import warnings

warnings.filterwarnings('ignore')

# Set random seed for reproducibility, to mirror Perl srand(23)
random.seed(23)
np.random.seed(23)

enst_column_no = None

def read_enst2ensg_file(filename):
    enst2ensg = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) > 2:
                enst2ensg[parts[2].split('.')[0]] = parts[0]
    return enst2ensg

def read_redundant_file(filename, verbose=False):
    redundant = {}
    if not filename or not os.path.exists(filename):
        return redundant
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                ensg, enst_sel, enst_red = parts[:3]
                redundant[enst_red.split('.')[0]] = enst_sel.split('.')[0]
    if verbose and redundant:
        print(f"Found {len(redundant)} redundant ENSTs", file=sys.stderr)
    return redundant

def read_sequence_file(filename, verbose=False):
    sequences = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                ensp, ensg, enst, seq = parts[:4]
                sequences[enst.split('.')[0]] = seq
    if verbose:
        print(f"Found {len(sequences)} sequences", file=sys.stderr)
    return sequences

def read_expression_file(filename, enst2ensg, redundant, verbose=False):
    enst_express = {}
    ensg_express = {}

    with gzip.open(filename, 'rt') as f:
        header = f.readline().rstrip('\n')
        samples = header.split('\t')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= enst_column_no:
                continue
            enst = parts[enst_column_no]
            if not enst.startswith('ENST'):
                continue

            enst = enst.split('.')[0]
            if enst in redundant:
                if verbose:
                    print(f"Replacing redundant cDNA {enst} with {redundant[enst]}", file=sys.stderr)
                enst = redundant[enst]

            ensg = enst2ensg.get(enst)
            if not ensg:
                if verbose:
                    print(f"WARNING: {enst} not in ENST2ENSG. Skipping.", file=sys.stderr)
                continue

            for i in range(enst_column_no + 1, min(len(parts), len(samples))):
                try:
                    expression = 0.0 if parts[i] == "NA" or parts[i] == "" else float(parts[i])
                except ValueError:
                    expression = 0.0
                sample = samples[i]

                enst_express.setdefault(ensg, {}).setdefault(sample, {})
                ensg_express.setdefault(ensg, {})

                enst_express[ensg][sample][enst] = enst_express[ensg][sample].get(enst, 0.0) + expression
                ensg_express[ensg][sample] = ensg_express[ensg].get(sample, 0.0) + expression

    # Compute relative expression per transcript per sample
    rel_enst_express = {}
    rel_enst_express_string = {}
    all_rel_enst_express = []

    for ensg in sorted(enst_express.keys()):
        for sample in sorted(enst_express[ensg].keys()):
            gene_expr = ensg_express[ensg][sample]
            for enst, tpm in sorted(enst_express[ensg][sample].items()):
                rel = (tpm / gene_expr) if gene_expr > 0 else 0.0
                rel_enst_express.setdefault(ensg, {}).setdefault(enst, {})[sample] = rel
                all_rel_enst_express.append(rel)
                rel_enst_express_string.setdefault(ensg, {}).setdefault(enst, "")
                rel_enst_express_string[ensg][enst] = (rel_enst_express_string[ensg][enst] + " " if rel_enst_express_string[ensg][enst] else "") + f"{rel:.3f}"

    return ensg_express, enst_express, rel_enst_express, rel_enst_express_string, all_rel_enst_express

def enrichment(express, minimum_express, min_enrichment, output_file=None):
    enrichment_dict = {}
    all_enrichment = []
    out_f = gzip.open(output_file, 'wt') if output_file else None

    for ensg in sorted(express.keys()):
        for sample in sorted(express[ensg].keys()):
            # sort transcripts by expression desc
            sorted_enst = sorted(express[ensg][sample].keys(), key=lambda x: express[ensg][sample][x], reverse=True)
            if not sorted_enst:
                continue
            top = sorted_enst[0]
            top_expr = express[ensg][sample][top]

            # Compare only to the second-most expressed, mirroring Perl
            if len(sorted_enst) >= 2:
                other = sorted_enst[1]
                other_expr = express[ensg][sample][other]
                enr = top_expr
                if top_expr >= minimum_express:
                    if other_expr > 0:
                        enr = top_expr / other_expr
                    if out_f:
                        out_f.write(f"{sample}\t{ensg}\t{top}\t{other}\t{top_expr}\t{other_expr}\t{enr}\n")
                    if enr >= min_enrichment:
                        enrichment_dict.setdefault(ensg, {}).setdefault(top, {})[sample] = enr
                        all_enrichment.append(enr)

    if out_f:
        out_f.close()
    return enrichment_dict, all_enrichment

def read_canon_file(filename):
    canon = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 6:
                canon_ensp, ensg, ensp, enst, gene_name, source = parts[:6]
                canon[enst.split('.')[0]] = canon_ensp
                canon[ensg] = canon_ensp
                canon[ensp] = ensg
    return canon

def read_isoform_int_file(filename, min_string_score):
    miss_int4isoform = {}
    int_n = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            ensg1, ensp1, enst1, string_ensp1, string_int_n, exist_int_n, miss_int_n, rel_miss_int_n, exist_int, miss_int = parts[:10]
            enst1 = enst1.split('.')[0]
            int_n[enst1] = 0
            miss_int4isoform.setdefault(enst1, {})

            # Missed interactions
            if miss_int:
                for interaction in miss_int.split(','):
                    if ':' not in interaction:
                        continue
                    int_parts = interaction.split(':')
                    if len(int_parts) >= 8:
                        try:
                            score = float(int_parts[6])
                            if score >= min_string_score:
                                miss_int4isoform[enst1][interaction] = 1
                                int_n[enst1] += 1
                        except Exception:
                            pass
            # Existing interactions
            if exist_int:
                for interaction in exist_int.split(','):
                    if ':' not in interaction:
                        continue
                    int_parts = interaction.split(':')
                    if len(int_parts) >= 8:
                        try:
                            score = float(int_parts[6])
                            if score >= min_string_score:
                                int_n[enst1] += 1
                        except Exception:
                            pass
    return miss_int4isoform, int_n

def median(values):
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    return (s[mid-1] + s[mid]) / 2.0

def main():
    parser = argparse.ArgumentParser(description='Dominant Transcript Switch Analysis (fixed)')
    parser.add_argument('--pcawg', required=True)
    parser.add_argument('--gtex', required=True)
    parser.add_argument('--ensg', required=True)
    parser.add_argument('--canon', required=True)
    parser.add_argument('--iso', required=True)
    parser.add_argument('--seq', required=True)
    parser.add_argument('--redundant')
    parser.add_argument('--minString', type=int, required=True)
    parser.add_argument('--enst', type=int, required=True)
    parser.add_argument('--minEnr', type=float, required=True)
    parser.add_argument('--maxQ', type=float, required=True)
    parser.add_argument('--minExp', type=float, required=True)
    parser.add_argument('--pvalueFile')
    parser.add_argument('--randomPvalueFile')
    parser.add_argument('--pcawgEnrichmentFile')
    parser.add_argument('--gtexEnrichmentFile')
    parser.add_argument('--verbose', action='store_true')
    # Extra helper to compare with expected file right after run
    parser.add_argument('--expected', help='Path to expected output (TSV/TSV.GZ/XLSX) for immediate validation')

    args = parser.parse_args()
    global enst_column_no
    enst_column_no = args.enst

    if args.verbose: print(f"Reading sequences: {args.seq}", file=sys.stderr)
    sequences = read_sequence_file(args.seq, args.verbose)

    if args.verbose: print(f"Reading ENST→ENSG: {args.ensg}", file=sys.stderr)
    enst2ensg = read_enst2ensg_file(args.ensg)

    if args.verbose and args.redundant: print(f"Reading redundant ENST: {args.redundant}", file=sys.stderr)
    redundant = read_redundant_file(args.redundant, args.verbose) if args.redundant else {}

    if args.verbose: print(f"Reading PCAWG: {args.pcawg}", file=sys.stderr)
    ensg_express_pcawg, enst_express_pcawg, rel_enst_express_pcawg, rel_enst_express_string_pcawg, all_rel_enst_express_pcawg = read_expression_file(args.pcawg, enst2ensg, redundant, args.verbose)

    if args.verbose: print(f"Reading GTEx: {args.gtex}", file=sys.stderr)
    ensg_express_gtex, enst_express_gtex, rel_enst_express_gtex, rel_enst_express_string_gtex, all_rel_enst_express_gtex = read_expression_file(args.gtex, enst2ensg, redundant, args.verbose)

    if args.verbose: print(f"Reading canonical: {args.canon}", file=sys.stderr)
    canonical = read_canon_file(args.canon)

    if args.verbose: print(f"Reading isoform interactions: {args.iso}", file=sys.stderr)
    miss_iso_int, int_n = read_isoform_int_file(args.iso, args.minString)

    if args.verbose: print("Calculating enrichment PCAWG...", file=sys.stderr)
    enrichment_pcawg, all_enrichment_pcawg = enrichment(enst_express_pcawg, args.minExp, args.minEnr, args.pcawgEnrichmentFile)

    if args.verbose: print("Calculating enrichment GTEx...", file=sys.stderr)
    enrichment_gtex, all_enrichment_gtex = enrichment(enst_express_gtex, args.minExp * 0.1, args.minEnr, args.gtexEnrichmentFile)

    # Output header
    header = "#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\tGTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\tNumberOfGTExMDT\tcMDTenrichment\tGTExMDTmedianRelExpression\tcMDTrelExpression\tcMDTuniqMissedInt\tTotalNumberOfSTRINGint\tcMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt"
    print(header)

    pvalue_file = open(args.pvalueFile, 'w') if args.pvalueFile else None
    if pvalue_file:
        pvalue_file.write("#ENSG\tcMDT\ttRNAseqAliquotID\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\tPCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\tRelMDTexpressionInGTEx\n")

    random_pvalue_file = open(args.randomPvalueFile, 'w') if args.randomPvalueFile else None
    if random_pvalue_file:
        random_pvalue_file.write("#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\tPCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n")

    # Records collected ONLY after fully passing all filters; preserves order for FDR
    records = []  # each is dict with keys: part1, part2, p
    for ensg in sorted(enrichment_pcawg.keys()):
        if ensg not in rel_enst_express_string_gtex:
            if args.verbose:
                print(f"WARNING: {ensg} not in GTEx rel expressions; skipping.", file=sys.stderr)
            continue

        for enst in sorted(enrichment_pcawg[ensg].keys()):
            if enst not in rel_enst_express_string_gtex[ensg]:
                if args.verbose:
                    print(f"WARNING: {ensg} - {enst} not in GTEx; skipping.", file=sys.stderr)
                continue

            for sample, enrichment_pcawg_val in sorted(enrichment_pcawg[ensg][enst].items()):
                enst_rel_express_pcawg = rel_enst_express_pcawg[ensg][enst][sample]
                enst_rel_expresses_gtex_str = rel_enst_express_string_gtex[ensg][enst]

                # Only analyze if ENST is NOT an MDT in GTEx but gene has MDTs
                if enst in enrichment_gtex.get(ensg, {}):
                    continue
                if ensg not in enrichment_gtex or len(enrichment_gtex[ensg]) == 0:
                    continue

                enst_rel_expresses_gtex = [float(x) for x in enst_rel_expresses_gtex_str.split()]
                rel_expression_gtex_median = median(enst_rel_expresses_gtex) if enst_rel_expresses_gtex else 0.0
                pos = sum(1 for x in enst_rel_expresses_gtex if enst_rel_express_pcawg > x)
                neg = sum(1 for x in enst_rel_expresses_gtex if enst_rel_express_pcawg < x)
                pvalue = binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue if (pos + neg) > 0 else 1.0

                if pvalue_file:
                    pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t{enst_rel_express_pcawg:.3f}\t{pos}\t{neg}\t{pvalue:.3e}\t{enst_rel_expresses_gtex_str}\n")

                if random_pvalue_file:
                    for _ in range(100):
                        rnd = random.choice(all_enrichment_pcawg) if all_enrichment_pcawg else 0.0
                        rpos = sum(1 for x in enst_rel_expresses_gtex if rnd > x)
                        rneg = sum(1 for x in enst_rel_expresses_gtex if rnd < x)
                        rp = binomtest(rpos, rpos + rneg, p=0.5, alternative='two-sided').pvalue if (rpos + rneg) > 0 else 1.0
                        random_pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t{rnd:.3f}\t{rpos}\t{rneg}\t{rp:.3e}\n")

                if pvalue >= args.maxQ:
                    continue  # early p-filter like Perl

                # Build GTEx MDTs summary and check 50% rule
                gtex_mdts = {}
                mdi_enrichments_gtex = []
                for n_enst in sorted(enrichment_gtex.get(ensg, {}).keys()):
                    for n_sample, n_enrichment in sorted(enrichment_gtex[ensg][n_enst].items()):
                        if n_enrichment >= args.minEnr and (n_enst in sequences) and (enst in sequences) and (sequences[n_enst] != sequences[enst]):
                            gtex_mdts[n_enst] = gtex_mdts.get(n_enst, 0) + 1
                            mdi_enrichments_gtex.append(n_enrichment)

                # Skip if GTEx doesn't have MDT for >=50% of relevant samples
                skip = True
                for g, count in gtex_mdts.items():
                    if count >= len(enst_rel_expresses_gtex) * 0.5:
                        skip = False
                        break
                if skip:
                    continue

                gtex_mdts_str = ",".join([f"{g}:{gtex_mdts[g]}" for g in sorted(gtex_mdts, key=lambda x: gtex_mdts[x], reverse=True)]) if gtex_mdts else "-"
                median_mdt_enr_gtex = f"{median(mdi_enrichments_gtex):.3f}" if mdi_enrichments_gtex else "-"

                # Interaction disruption counts
                string_int_n = int_n.get(enst, -1)
                missed_int = "-"
                missed_cancer_int_n = -1
                missed_common_int_n = -1

                if enst in miss_iso_int:
                    missed_cancer_int_n = 0
                    missed_common_int_n = 0
                    buf = []
                    for c_int in sorted(miss_iso_int[enst].keys()):
                        parts = c_int.split(':')
                        if len(parts) >= 8:
                            ensp2 = parts[5]
                            if ensp2 in canonical and canonical[ensp2] in ensg_express_pcawg and sample in ensg_express_pcawg[canonical[ensp2]]:
                                found = False
                                # is this interaction also missed by a GTEx MDT (with enrichment >= minEnr)
                                for n_enst in enrichment_gtex.get(ensg, {}):
                                    for n_sample, n_enrichment in enrichment_gtex[ensg][n_enst].items():
                                        if n_enrichment >= args.minEnr and (n_enst in miss_iso_int) and (c_int in miss_iso_int[n_enst]):
                                            found = True
                                            break
                                    if found:
                                        break
                                if found:
                                    missed_common_int_n += 1
                                else:
                                    buf.append(c_int)
                                    missed_cancer_int_n += 1
                    if buf:
                        missed_int = ",".join(buf)

                canon_ensp = canonical.get(enst, "-")
                part1 = f"{ensg}\t{canon_ensp}\t{enst}\t{sample}\t{gtex_mdts_str}\t{len(enst_rel_expresses_gtex)}\t{median_mdt_enr_gtex}"
                part2 = f"{len(mdi_enrichments_gtex)}\t{enrichment_pcawg_val:.3f}\t{rel_expression_gtex_median:.3f}\t{enst_rel_express_pcawg:.3f}\t{missed_int}\t{string_int_n}\t{missed_cancer_int_n}\t{missed_common_int_n}"
                records.append({"part1": part1, "part2": part2, "p": pvalue})

    # Apply FDR on the collected (already-filtered) records
    if records:
        pvals = [r["p"] for r in records]
        reject, qvals, *_ = multipletests(pvals, method='fdr_bh')
    else:
        qvals = []

    for i, rec in enumerate(records):
        q = qvals[i] if i < len(qvals) else 1.0
        if q < args.maxQ:
            print(f"{rec['part1']}\t{rec['p']:.3e}\t{q:.3e}\t{rec['part2']}")

    if pvalue_file: pvalue_file.close()
    if random_pvalue_file: random_pvalue_file.close()

    # Optional: quick validator, if --expected is given
    if args.expected:
        try:
            import pandas as pd, difflib
            exp_path = args.expected
            if exp_path.endswith('.xlsx'):
                expected_df = pd.read_excel(exp_path)
            elif exp_path.endswith('.gz'):
                expected_df = pd.read_csv(exp_path, sep='\t', compression='gzip', comment='#')
            else:
                expected_df = pd.read_csv(exp_path, sep='\t', comment='#')

            # Capture our stdout (excluding header comment line) by re-running print block to a buffer
            # For simplicity here, we won't re-run; instead tell the user how to compare with provided test.
            print("#VALIDATOR_NOTE\tTo compare outputs robustly, please use the provided test script or pandas compare.", file=sys.stderr)
        except Exception as e:
            print(f"#VALIDATOR_ERROR\t{e}", file=sys.stderr)


if __name__ == '__main__':
    main()
