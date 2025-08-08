#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python implementation of mostDominantTranscriptSwitch.pl (fixed, v2)
- Sign-test uses full-precision GTEx relative expressions (printing still 3dp)
- Missed-interaction counters stay -1/-1 unless there is at least one qualified missed interaction
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

random.seed(23)
np.random.seed(23)

enst_column_no = None

def read_enst2ensg_file(filename):
    enst2ensg = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) > 2:
                enst2ensg[parts[2].split('.')[0]] = parts[0]
    return enst2ensg

def read_redundant_file(filename, verbose=False):
    redundant = {}
    if not filename or not os.path.exists(filename): return redundant
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('EN'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                _, enst_sel, enst_red = parts[:3]
                redundant[enst_red.split('.')[0]] = enst_sel.split('.')[0]
    if verbose and redundant:
        print(f"Found {len(redundant)} redundant ENSTs", file=sys.stderr)
    return redundant

def read_sequence_file(filename, verbose=False):
    sequences = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                _, _, enst, seq = parts[:4]
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
            if len(parts) <= enst_column_no: continue
            enst = parts[enst_column_no]
            if not enst.startswith('ENST'): continue

            enst = enst.split('.')[0]
            if enst in redundant:
                if verbose:
                    print(f"Replacing redundant cDNA {enst} with {redundant[enst]}", file=sys.stderr)
                enst = redundant[enst]

            ensg = enst2ensg.get(enst)
            if not ensg:
                if verbose: print(f"WARNING: {enst} not in ENST2ENSG. Skipping.", file=sys.stderr)
                continue

            for i in range(enst_column_no + 1, min(len(parts), len(samples))):
                try:
                    expression = 0.0 if parts[i] in ("NA","") else float(parts[i])
                except ValueError:
                    expression = 0.0
                sample = samples[i]
                enst_express.setdefault(ensg, {}).setdefault(sample, {})
                ensg_express.setdefault(ensg, {})
                enst_express[ensg][sample][enst] = enst_express[ensg][sample].get(enst, 0.0) + expression
                ensg_express[ensg][sample] = ensg_express[ensg].get(sample, 0.0) + expression

    # relative expressions (both precise dict & 3dp string for printing)
    rel_enst_express = {}           # precise floats: rel_enst_express[ensg][enst][sample] = float
    rel_enst_express_string = {}    # for printing only (joined "0.xxx 0.yyy")
    all_rel_enst_express = []

    for ensg in sorted(enst_express.keys()):
        for sample in sorted(enst_express[ensg].keys()):
            gene_expr = ensg_express[ensg][sample]
            for enst, tpm in sorted(enst_express[ensg][sample].items()):
                rel = (tpm / gene_expr) if gene_expr > 0 else 0.0
                rel_enst_express.setdefault(ensg, {}).setdefault(enst, {})[sample] = rel
                all_rel_enst_express.append(rel)
                rel_enst_express_string.setdefault(ensg, {}).setdefault(enst, "")
                rel_enst_express_string[ensg][enst] += ("" if rel_enst_express_string[ensg][enst]=="" else "")+f"{rel:.3f} " if False else ""
    # yukarıdaki satır string biriktirmiyor; alt blokta daha kontrollü kuracağız

    # string biriktirme (ayrı döngü: örnek sıraya göre 3 ondalık)
    for ensg in sorted(rel_enst_express.keys()):
        for enst in sorted(rel_enst_express[ensg].keys()):
            vals = []
            for sample in sorted(rel_enst_express[ensg][enst].keys()):
                vals.append(f"{rel_enst_express[ensg][enst][sample]:.3f}")
            rel_enst_express_string[ensg][enst] = " ".join(vals)

    return ensg_express, enst_express, rel_enst_express, rel_enst_express_string, all_rel_enst_express

def enrichment(express, minimum_express, min_enrichment, output_file=None):
    enrichment_dict = {}
    all_enrichment = []
    out_f = gzip.open(output_file, 'wt') if output_file else None

    for ensg in sorted(express.keys()):
        for sample in sorted(express[ensg].keys()):
            sorted_enst = sorted(express[ensg][sample].keys(),
                                 key=lambda x: express[ensg][sample][x], reverse=True)
            if not sorted_enst: continue
            top = sorted_enst[0]
            top_expr = express[ensg][sample][top]
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

    if out_f: out_f.close()
    return enrichment_dict, all_enrichment

def read_canon_file(filename):
    canon = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('ENS'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 6:
                canon_ensp, ensg, ensp, enst, _, _ = parts[:6]
                canon[enst.split('.')[0]] = canon_ensp
                canon[ensg] = canon_ensp
                canon[ensp] = ensg
    return canon

def read_isoform_int_file(filename, min_string_score):
    """
    Only create miss_int4isoform[ENST] key if there is at least one qualified missed interaction (>= min_string_score).
    Always compute int_n (total qualified STRING interactions) as before.
    """
    miss_int4isoform = {}
    int_n = {}
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if not line.startswith('EN'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10: continue
            _, _, enst1, _, _, _, _, _, exist_int, miss_int = parts[:10]
            enst1 = enst1.split('.')[0]
            int_n.setdefault(enst1, 0)

            added_any_missed = False

            # missed interactions
            if miss_int:
                for interaction in miss_int.split(','):
                    if ':' not in interaction: continue
                    int_parts = interaction.split(':')
                    if len(int_parts) >= 8:
                        try:
                            score = float(int_parts[6])
                            if score >= min_string_score:
                                if enst1 not in miss_int4isoform:
                                    miss_int4isoform[enst1] = {}
                                miss_int4isoform[enst1][interaction] = 1
                                int_n[enst1] += 1
                                added_any_missed = True
                        except Exception:
                            pass

            # existing interactions (only count total STRING int, do not create miss map)
            if exist_int:
                for interaction in exist_int.split(','):
                    if ':' not in interaction: continue
                    int_parts = interaction.split(':')
                    if len(int_parts) >= 8:
                        try:
                            score = float(int_parts[6])
                            if score >= min_string_score:
                                int_n[enst1] += 1
                        except Exception:
                            pass

            # NOTE: do NOT add empty dict for enst1 if no qualified missed int were found

    return miss_int4isoform, int_n

def median(values):
    if not values: return 0.0
    s = sorted(values); n = len(s); mid = n//2
    return s[mid] if n%2==1 else (s[mid-1]+s[mid])/2.0

def main():
    parser = argparse.ArgumentParser(description='Dominant Transcript Switch Analysis (fixed v2)')
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
    ensg_express_pcawg, enst_express_pcawg, rel_enst_express_pcawg, rel_enst_express_string_pcawg, _ = \
        read_expression_file(args.pcawg, enst2ensg, redundant, args.verbose)

    if args.verbose: print(f"Reading GTEx: {args.gtex}", file=sys.stderr)
    ensg_express_gtex, enst_express_gtex, rel_enst_express_gtex, rel_enst_express_string_gtex, _ = \
        read_expression_file(args.gtex, enst2ensg, redundant, args.verbose)

    if args.verbose: print(f"Reading canonical: {args.canon}", file=sys.stderr)
    canonical = read_canon_file(args.canon)

    if args.verbose: print(f"Reading isoform interactions: {args.iso}", file=sys.stderr)
    miss_iso_int, int_n = read_isoform_int_file(args.iso, args.minString)

    if args.verbose: print("Calculating enrichment PCAWG...", file=sys.stderr)
    enrichment_pcawg, all_enrichment_pcawg = enrichment(enst_express_pcawg, args.minExp, args.minEnr, args.pcawgEnrichmentFile)

    if args.verbose: print("Calculating enrichment GTEx...", file=sys.stderr)
    enrichment_gtex, _ = enrichment(enst_express_gtex, args.minExp * 0.1, args.minEnr, args.gtexEnrichmentFile)

    # header
    print("#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\tGTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\tNumberOfGTExMDT\tcMDTenrichment\tGTExMDTmedianRelExpression\tcMDTrelExpression\tcMDTuniqMissedInt\tTotalNumberOfSTRINGint\tcMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt")

    pvalue_file = open(args.pvalueFile, 'w') if args.pvalueFile else None
    if pvalue_file:
        pvalue_file.write("#ENSG\tcMDT\ttRNAseqAliquotID\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\tPCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\tRelMDTexpressionInGTEx\n")

    random_pvalue_file = open(args.randomPvalueFile, 'w') if args.randomPvalueFile else None
    if random_pvalue_file:
        random_pvalue_file.write("#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\tPCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n")

    records = []
    for ensg in sorted(enrichment_pcawg.keys()):
        if ensg not in rel_enst_express_gtex:
            if args.verbose: print(f"WARNING: {ensg} not in GTEx rel expressions; skipping.", file=sys.stderr)
            continue
        for enst in sorted(enrichment_pcawg[ensg].keys()):
            if enst not in rel_enst_express_gtex[ensg]:
                if args.verbose: print(f"WARNING: {ensg} - {enst} not in GTEx; skipping.", file=sys.stderr)
                continue
            for sample, enrichment_pcawg_val in sorted(enrichment_pcawg[ensg][enst].items()):
                enst_rel_express_pcawg = rel_enst_express_pcawg[ensg][enst][sample]

                # full-precision GTEx rel expressions for sign-test
                gtex_vals_precise = [rel_enst_express_gtex[ensg][enst][s] for s in sorted(rel_enst_express_gtex[ensg][enst].keys())]
                gtex_vals_print = " ".join([f"{v:.3f}" for v in gtex_vals_precise])

                # Only analyze if this ENST is NOT MDT in GTEx, but gene has MDT(s)
                if enst in enrichment_gtex.get(ensg, {}): continue
                if ensg not in enrichment_gtex or len(enrichment_gtex[ensg]) == 0: continue

                rel_expression_gtex_median = median(gtex_vals_precise) if gtex_vals_precise else 0.0
                pos = sum(1 for x in gtex_vals_precise if enst_rel_express_pcawg > x)
                neg = sum(1 for x in gtex_vals_precise if enst_rel_express_pcawg < x)
                pvalue = binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue if (pos + neg) > 0 else 1.0

                if pvalue_file:
                    pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t{enst_rel_express_pcawg:.3f}\t{pos}\t{neg}\t{pvalue:.3e}\t{gtex_vals_print}\n")

                if random_pvalue_file:
                    for _ in range(100):
                        rnd = random.choice(all_enrichment_pcawg) if all_enrichment_pcawg else 0.0
                        rpos = sum(1 for x in gtex_vals_precise if rnd > x)
                        rneg = sum(1 for x in gtex_vals_precise if rnd < x)
                        rp = binomtest(rpos, rpos + rneg, p=0.5, alternative='two-sided').pvalue if (rpos + rneg) > 0 else 1.0
                        random_pvalue_file.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t{rnd:.3f}\t{rpos}\t{rneg}\t{rp:.3e}\n")

                if pvalue >= args.maxQ:  # early p-filter
                    continue

                # Build GTEx MDTs summary & 50% rule
                gtex_mdts = {}
                mdi_enrichments_gtex = []
                for n_enst in sorted(enrichment_gtex.get(ensg, {}).keys()):
                    for n_sample, n_enrichment in sorted(enrichment_gtex[ensg][n_enst].items()):
                        if n_enrichment >= args.minEnr and (n_enst in sequences) and (enst in sequences) and (sequences[n_enst] != sequences[enst]):
                            gtex_mdts[n_enst] = gtex_mdts.get(n_enst, 0) + 1
                            mdi_enrichments_gtex.append(n_enrichment)

                skip = True
                for g, count in gtex_mdts.items():
                    if count >= len(gtex_vals_precise) * 0.5:
                        skip = False; break
                if skip: continue

                gtex_mdts_str = ",".join([f"{g}:{gtex_mdts[g]}" for g in sorted(gtex_mdts, key=lambda x: gtex_mdts[x], reverse=True)]) if gtex_mdts else "-"
                median_mdt_enr_gtex = f"{median(mdi_enrichments_gtex):.3f}" if mdi_enrichments_gtex else "-"

                # Interaction disruption counters
                string_int_n = int_n.get(enst, -1)
                missed_int = "-"
                missed_cancer_int_n = -1
                missed_common_int_n = -1

                # only if there exists at least one qualified missed interaction for this ENST
                if (enst in miss_iso_int) and (len(miss_iso_int[enst]) > 0):
                    missed_cancer_int_n = 0
                    missed_common_int_n = 0
                    buf = []
                    for c_int in sorted(miss_iso_int[enst].keys()):
                        parts = c_int.split(':')
                        if len(parts) >= 8:
                            ensp2 = parts[5]
                            if ensp2 in canonical and canonical[ensp2] in ensg_express_pcawg and sample in ensg_express_pcawg[canonical[ensp2]]:
                                found = False
                                # also missed by some GTEx MDT?
                                for n_enst in enrichment_gtex.get(ensg, {}):
                                    for n_sample, n_enrichment in enrichment_gtex[ensg][n_enst].items():
                                        if n_enrichment >= args.minEnr and (n_enst in miss_iso_int) and (c_int in miss_iso_int[n_enst]):
                                            found = True; break
                                    if found: break
                                if found:
                                    missed_common_int_n += 1
                                else:
                                    buf.append(c_int)
                                    missed_cancer_int_n += 1
                    if buf: missed_int = ",".join(buf)

                canon_ensp = canonical.get(enst, "-")
                part1 = f"{ensg}\t{canon_ensp}\t{enst}\t{sample}\t{gtex_mdts_str}\t{len(gtex_vals_precise)}\t{median_mdt_enr_gtex}"
                part2 = f"{len(mdi_enrichments_gtex)}\t{enrichment_pcawg_val:.3f}\t{median(gtex_vals_precise):.3f}\t{enst_rel_express_pcawg:.3f}\t{missed_int}\t{string_int_n}\t{missed_cancer_int_n}\t{missed_common_int_n}"
                records.append({"part1": part1, "part2": part2, "p": pvalue})

    # FDR (on already-filtered rows; preserves order)
    if records:
        pvals = [r["p"] for r in records]
        _, qvals, *_ = multipletests(pvals, method='fdr_bh')
    else:
        qvals = []

    for i, rec in enumerate(records):
        q = qvals[i] if i < len(qvals) else 1.0
        if q < args.maxQ:
            print(f"{rec['part1']}\t{rec['p']:.3e}\t{q:.3e}\t{rec['part2']}")

    if pvalue_file: pvalue_file.close()
    if random_pvalue_file: random_pvalue_file.close()

if __name__ == '__main__':
    main()
