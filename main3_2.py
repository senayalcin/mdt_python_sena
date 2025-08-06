#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Abdullah Kahraman (translated+annotated)
Date:   06.02.2018 (original Perl)
Purpose: Determines cancer-specific MDT and compares with GTEx using Sign-test.
"""

import argparse, sys, gzip, os, statistics
from collections import defaultdict
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def parse_args():
    p = argparse.ArgumentParser(
        description="Determines cancer-specific MDT based on 2-fold expression difference "
                    "to minor dominant transcripts and compares with GTEx using Sign-test."
    )
    p.add_argument('-pcawg',        '--pcawg',        required=True,  help='PCAWG transcript TPM file (gzipped)')
    p.add_argument('-gtex',         '--gtex',         required=True,  help='GTEx transcript TPM file (gzipped)')
    p.add_argument('-ensg',         '--ensg',         required=True,  help='ENST→ENSG mapping file (gzipped)')
    p.add_argument('-canon',        '--canon',        required=True,  help='Canonical isoforms file')
    p.add_argument('-iso',          '--iso',          required=True,  help='Isoform interaction file')
    p.add_argument('-seq',          '--seq',          required=True,  help='ENSP sequence file')
    p.add_argument('-redundant',    '--redundant',    required=True,  help='Redundant ENST replacement file')
    p.add_argument('-minString',    '--minString',    type=int,    default=900,   help='Minimum STRING score')
    p.add_argument('-enstColumnNo', '--enstColumnNo', type=int,    default=0,     help='ENST column index (0-based)')
    p.add_argument('-minEnr',       '--minEnr',       type=float,  default=2.0,   help='Minimum enrichment')
    p.add_argument('-maxQ',         '--maxQ',         type=float,  default=0.01,  help='Maximum FDR q-value')
    p.add_argument('-minExp',       '--minExp',       type=float,  default=2.0,   help='Minimum expression')
    p.add_argument('-v',            '--verbose',      action='store_true',        help='Verbose progress')
    return p.parse_args()

args = parse_args()



def open_maybe_gzip(path, mode='rt'):
    return gzip.open(path, mode=mode) if path.endswith('.gz') else open(path, mode=mode, encoding='utf-8', errors='ignore')



def read_enst2ensg_file(path):
    mapping = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('ENS'): continue
            f = line.rstrip('\n').split('\t')
            mapping[f[2]] = f[0]
    if args.verbose:
        print(f"Loaded {len(mapping)} ENST→ENSG mappings", file=sys.stderr)
    return mapping


def read_redundant_file(path):
    redundant = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'): continue
            p = line.rstrip('\n').split('\t')
            redundant[p[2]] = p[1]
    if args.verbose:
        print(f"Loaded {len(redundant)} redundant transcript mappings", file=sys.stderr)
    return redundant


def read_canon_file(path):
    canon = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('ENS'): continue
            p = line.rstrip('\n').split('\t')
            # ENST→ENSPcanon, ENSG→ENSPcanon, ENSP→ENSG
            canon[p[3]] = p[0]
            canon[p[1]] = p[0]
            canon[p[2]] = p[1]
    if args.verbose:
        print(f"Loaded {len(canon)} canonical entries", file=sys.stderr)
    return canon


def read_isoform_int_file(path, min_string_score):
    miss_int= {}
    exist_counts = {}
    int_counts = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'): continue
            parts = line.rstrip('\n').split('\t')
            raw_id = parts[2]
            enst1  = raw_id.split('.')[0]      # versiyonsuz
            raw    = parts[-1]                 # ← MISSINTS sütunu (eksik etkileşimler)
            exist_counts[enst1] = int(parts[5])  # <--- pull in ExistIntN directly
            int_counts[enst1] = 0
            for entry in raw.split(','):
                if ':' not in entry: continue
                fields = entry.split(':')
                score  = float(fields[-2])
                if score < min_string_score: continue
                miss_int.setdefault(enst1, {})[entry] = True
                int_counts[enst1] += 1
    if args.verbose:
        print(f"Processed interactions: {len(miss_int)} transcripts with STRING≥{args.minString}", file=sys.stderr)
    return miss_int, int_counts, exist_counts


def read_sequence_file(path):
    sequences = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'): continue
            p = line.rstrip('\n').split('\t')
            sequences[p[2]] = p[3]
    if args.verbose:
        print(f"Loaded {len(sequences)} transcript sequences", file=sys.stderr)
    return sequences


def read_expression_file(path):
    enst2ensg = read_enst2ensg_file(args.ensg)
    redundant = read_redundant_file(args.redundant)
    enst_expr = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    ensg_expr = defaultdict(lambda: defaultdict(float))
    rel_expr, rel_str, all_rel = defaultdict(lambda: defaultdict(dict)), defaultdict(dict), []
    with open_maybe_gzip(path) as fh:
        header = next(fh).rstrip('\n').split('\t')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            enst = parts[args.enstColumnNo]
            if enst in redundant:
                if args.verbose:
                    print(f"Replacing redundant {enst} → {redundant[enst]}", file=sys.stderr)
                enst = redundant[enst]
            if not enst.startswith('ENST0'): continue
            enst = enst.split('.')[0]
            if enst not in enst2ensg:
                print(f"WARNING: {enst} not in ENST→ENSG mapping", file=sys.stderr)
                continue
            ensg = enst2ensg[enst]
            for i in range(args.enstColumnNo+1, len(parts)):
                sample, val = header[i], parts[i]
                expr = 0.0 if val == 'NA' else float(val)
                enst_expr[ensg][sample][enst] += expr
                ensg_expr[ensg][sample]   += expr
    for gene, samp_dict in enst_expr.items():
        for sample, iso_dict in samp_dict.items():
            total = ensg_expr[gene][sample]
            for iso, val in iso_dict.items():
                rel = val/total if total>0 else 0.0
                rel_expr[gene][iso][sample] = rel
                all_rel.append(rel)
                prev = rel_str[gene].get(iso, "")
                rel_str[gene][iso] = (prev+" "+f"{rel:.3f}") if prev else f"{rel:.3f}"
    if args.verbose:
        print(f"Processed expression: {len(ensg_expr)} genes, {len(all_rel)} isoform rel-values", file=sys.stderr)
    return ensg_expr, enst_expr, rel_expr


def enrichment(data, min_exp):
    enrich, all_e = {}, []
    for gene, samp_dict in data.items():
        for sample, iso_dict in samp_dict.items():
            sorted_iso = sorted(iso_dict.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_iso) < 2: continue
            top_iso, top_val = sorted_iso[0]
            second_iso, sec_val = sorted_iso[1]
            if top_val < min_exp: continue
            enr = top_val/sec_val if sec_val>0 else float('inf')
            if enr >= args.minEnr:
                enrich.setdefault(gene, {}).setdefault(top_iso, {})[sample] = enr
                all_e.append(enr)
    if args.verbose:
        print(f"Computed enrichment: {len(all_e)} events across {len(enrich)} genes", file=sys.stderr)
    return enrich


def main():
    sequences            = read_sequence_file(args.seq)
    miss_int, int_counts, exist_counts = read_isoform_int_file(args.iso, args.minString)
    canon_map            = read_canon_file(args.canon)

    ensg_pc, enst_pc, rel_pc = read_expression_file(args.pcawg)
    ensg_gt, enst_gt, rel_gt = read_expression_file(args.gtex)

    enrich_pc = enrichment(enst_pc, args.minExp)
    enrich_gt = enrichment(enst_gt, args.minExp)

    output1, output2, pvals = [], [], []

    for gene in sorted(enrich_pc):
        if gene not in rel_gt: continue
        if gene not in enrich_gt or not enrich_gt[gene]: continue

        for iso in sorted(enrich_pc[gene]):
            if iso not in rel_gt[gene]: continue
            if iso in enrich_gt[gene]: continue

            # ← Burası değişti: sorted örnekler
            for sample in sorted(enrich_pc[gene][iso]):
                enr_pc = enrich_pc[gene][iso][sample]
                vals_gt = list(rel_gt[gene][iso].values())
                if not vals_gt: continue

                med_gt  = statistics.median(vals_gt)
                rel_c   = rel_pc[gene][iso][sample]
                pos     = sum(1 for v in vals_gt if rel_c > v)
                neg     = sum(1 for v in vals_gt if rel_c < v)
                if pos + neg == 0:
                    continue
                pv = binomtest(pos, pos+neg, p=0.5, alternative='two-sided').pvalue

                # GTEx MDT analysis
                gtex_counts, mdi_vals = {}, []
                for niso in enrich_gt[gene]:
                    for ns in enrich_gt[gene][niso]:
                        n_enr = enrich_gt[gene][niso][ns]
                        if n_enr >= args.minEnr and sequences.get(niso) and sequences[niso] != sequences[iso]:
                            gtex_counts[niso] = gtex_counts.get(niso, 0) + 1
                            mdi_vals.append(n_enr)

                total_gt = len(vals_gt)
                skip     = True
                gtex_list = []
                for niso, cnt in sorted(gtex_counts.items(), key=lambda x: x[1], reverse=True):
                    if cnt >= total_gt * 0.5:
                        skip = False
                    gtex_list.append(f"{niso}:{cnt}")
                if skip:
                    continue

                gtexMDTs   = ",".join(gtex_list) if gtex_list else "-"
                med_mdi_gt = statistics.median(mdi_vals) if mdi_vals else "-"

                # Interaction disruption
                # string_n = -1
                # if iso in int_counts:
                #     string_n = int_counts[iso]
                string_n = exist_counts.get(iso, -1)

                missed_can = missed_com = -1
                missed_list = []
                if iso in miss_int:
                    missed_can = missed_com = 0
                for c in miss_int.get(iso, {}):
                    parts = c.split(':')
                    if len(parts) < 7: continue
                    ensp2 = parts[5]
                    g2 = canon_map.get(ensp2)
                    if g2 in ensg_pc and sample in ensg_pc[g2]:
                        found = any(c in miss_int.get(niso, {}) for niso in enrich_gt[gene])
                        if found:
                            missed_com += 1
                        else:
                            missed_can += 1; missed_list.append(c)
                missed_str = ",".join(missed_list) if missed_list else "-"
                can_ensp   = canon_map.get(iso, "-")
                num_gt     = len(vals_gt)

                o1 = (
                    f"{gene}\t{can_ensp}\t{iso}\t{sample}\t"
                    f"{gtexMDTs}\t{num_gt}\t{med_mdi_gt:.3f}"
                )
                o2 = f"{len(mdi_vals)}\t{enr_pc:.3f}\t{med_gt:.3f}\t{rel_c:.3f}\t{missed_str}\t{string_n}\t{missed_can}\t{missed_com}"
                output1.append(o1)
                output2.append(o2)
                pvals.append(pv)

    if not pvals:
        if args.verbose:
            print("No P-values to adjust; exiting.", file=sys.stderr)
        sys.exit(0)

    reject, qvals, _, _ = multipletests(pvals, alpha=args.maxQ, method='fdr_bh')

    header = [
        "#ENSG",
        "ENSPcanon",
        "cMDT",
        "RNAseqAliquotID",
        "GTExMDT",
        "TotalRelevantGTExSamples",
        "GTExMDTmedianEnrichment",
        "Pvalue",
        "Qvalue",
        "NumberOfGTExMDT",
        "cMDTenrichment",
        "GTExMDTmedianRelExpression",
        "cMDTrelExpression",
        "cMDTuniqMissedInt",
        "TotalNumberOfSTRINGint",
        "cMDTnumberOfUniqMissedInt",
        "NumberOfCommonMissedInt"
    ]
    return header, output1, pvals, qvals, output2


if __name__ == '__main__':
    header, rows, pvals, qvals, output2 = main()

    out_path = 'output_py/python_main3_2.tsv.gz'
    with gzip.open(out_path, 'wt') as out_fh:
        out_fh.write("\t".join(header) + "\n")
        for o1, pv, qv, o2 in zip(rows, pvals, qvals, output2):
            out_fh.write(f"{o1}\t{pv:.3e}\t{qv:.3e}\t{o2}\n")
    if args.verbose:
        print(f"Written: {out_path}", file=sys.stderr)