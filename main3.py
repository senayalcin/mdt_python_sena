#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Abdullah Kahraman (translated to Python)
Date:   06.02.2018 (original Perl)

Purpose:
    Determines cancer-specific MDT based on 2-fold expression difference to
    minor dominant transcripts and compares the relative expression difference
    to expression differences in GTEx using a Sign-test. SENA YALÇINNNNNNNN

Usage example (PowerShell, single line):
    python3 main3.py \
      -pcawg Datas/pcawg.rnaseq.transcript.expr.tpm.tsv.gz \
      -gtex  Datas/GTEX_v4.pcawg.transcripts.tpm.tsv.gz \
      -ensg  Datas/ensg_ensp_enst_ense_geneName_v75.tsv.gz \
      -canon Datas/canonEnsp_ensg_ensp_enst_geneName_v75.tsv.gz \
      -iso   Datas/interactionsInIsoforms_900_2.tsv.gz \
      -seq   Datas/ensp_ensg_enst_sequence.tsv.gz \
      -redundant Datas/redundantENST_v75.txt \
      -minString 900 -minEnr 2 -maxQ 0.01 -minExp 2 -enstColumnNo 1 \
      -v \
      -pcawgEnrichmentFile pcawg_enrich.tsv \
      -gtexEnrichmentFile  gtex_enrich.tsv
"""
import argparse
import sys
import gzip
import statistics
from collections import defaultdict

from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Determines cancer-specific MDT based on 2-fold expression difference "
            "to minor dominant transcripts and compares with GTEx using Sign-test."
        )
    )
    parser.add_argument('-pcawg', '--pcawg',    required=True,
                        help='PCAWG transcript TPM file (gzipped)')
    parser.add_argument('-gtex',  '--gtex',     required=True,
                        help='GTEx transcript TPM file (gzipped)')
    parser.add_argument('-ensg',  '--ensg',     required=True,
                        help='ENST→ENSG mapping file (gzipped)')
    parser.add_argument('-canon', '--canon',    required=True,
                        help='Canonical isoforms file (plain or gzipped)')
    parser.add_argument('-iso',   '--iso',      required=True,
                        help='Isoform interaction file (plain or gzipped)')
    parser.add_argument('-seq',   '--seq',      required=True,
                        help='ENSP sequence file (plain or gzipped)')
    parser.add_argument('-redundant', '--redundant', required=True,
                        help='Redundant ENST replacement file')
    parser.add_argument('-minString', '--minString', type=int, default=900,
                        help='Minimum STRING score [default: 900]')
    parser.add_argument('-enstColumnNo', '--enstColumnNo', type=int, default=0,
                        help='ENST column index in expression files (0-based)')
    parser.add_argument('-minEnr', '--minEnr',   type=float, default=2.0,
                        help='Minimum enrichment [default: 2.0]')
    parser.add_argument('-maxQ',   '--maxQ',     type=float, default=0.01,
                        help='Maximum FDR q-value [default: 0.01]')
    parser.add_argument('-minExp', '--minExp',   type=float, default=2.0,
                        help='Minimum expression [default: 2.0]')
    parser.add_argument('-v',      '--verbose', action='store_true',
                        help='Print progress messages to stderr')
    parser.add_argument('-pcawgEnrichmentFile', '--pcawgEnrichmentFile',
                        help='PCAWG enrichment output file', default=None)
    parser.add_argument('-gtexEnrichmentFile',  '--gtexEnrichmentFile',
                        help='GTEx enrichment output file',  default=None)
    return parser.parse_args()


args = parse_args()


def open_maybe_gzip(path, mode='rt'):
    if path.endswith('.gz'):
        return gzip.open(path, mode=mode)
    return open(path, mode=mode, encoding='utf-8', errors='ignore')


def read_enst2ensg_file(path):
    mapping = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('ENS'):
                continue
            fields = line.rstrip('\n').split('\t')
            mapping[fields[2]] = fields[0]
    if args.verbose:
        print(f"Loaded {len(mapping)} ENST→ENSG mappings", file=sys.stderr)
    return mapping


def read_redundant_file(path):
    redundant = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            redundant[parts[2]] = parts[1]
    if args.verbose:
        print(f"Loaded {len(redundant)} redundant transcript mappings", file=sys.stderr)
    return redundant


def read_canon_file(path):
    canon = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('ENS'):
                continue
            parts = line.rstrip('\n').split('\t')
            canon_ensp, ensg, ensp, enst = parts[0], parts[1], parts[2], parts[3]
            canon[enst] = canon_ensp
            canon[ensg] = canon_ensp
            canon[ensp] = ensg
    if args.verbose:
        print(f"Loaded {len(canon)} canonical entries", file=sys.stderr)
    return canon


def read_isoform_int_file(path, min_string_score):
    miss_int = {}
    int_counts = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            enst1 = parts[2]
            raw_miss = parts[-1]
            int_counts[enst1] = 0
            for entry in raw_miss.split(','):
                if ':' not in entry:
                    continue
                fields = entry.split(':')
                score = float(fields[-2])
                if score < min_string_score:
                    continue
                miss_int.setdefault(enst1, {})[entry] = True
                int_counts[enst1] += 1
    if args.verbose:
        print(f"Processed interactions: {len(miss_int)} transcripts with STRING≥{min_string_score}", file=sys.stderr)
    return miss_int, int_counts


def read_sequence_file(path):
    sequences = {}
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line.startswith('EN'):
                continue
            parts = line.rstrip('\n').split('\t')
            sequences[parts[2]] = parts[3]
    if args.verbose:
        print(f"Loaded {len(sequences)} transcript sequences", file=sys.stderr)
    return sequences


def read_expression_file(path, sequences):
    enst2ensg = read_enst2ensg_file(args.ensg)
    redundant = read_redundant_file(args.redundant)

    # nested structures with defaultdict
    enst_expr = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    ensg_expr = defaultdict(lambda: defaultdict(float))
    rel_expr  = defaultdict(lambda: defaultdict(dict))
    rel_str   = defaultdict(dict)
    all_rel   = []

    with open_maybe_gzip(path) as fh:
        header = next(fh).rstrip('\n').split('\t')
        samples = header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            enst = parts[args.enstColumnNo]
            if enst in redundant:
                if args.verbose:
                    print(f"Replacing redundant {enst} with {redundant[enst]}", file=sys.stderr)
                enst = redundant[enst]
            if not enst.startswith('ENST0'):
                continue
            enst = enst.split('.')[0]
            if enst not in enst2ensg:
                print(f"WARNING: {enst} not in ENST→ENSG mapping", file=sys.stderr)
                continue
            ensg = enst2ensg[enst]
            for i in range(args.enstColumnNo + 1, len(parts)):
                sample = samples[i]
                val    = parts[i]
                expr   = 0.0 if val == 'NA' else float(val)
                enst_expr[ensg][sample][enst] += expr
                ensg_expr[ensg][sample]     += expr

    for gene, samp_dict in enst_expr.items():
        for sample, iso_dict in samp_dict.items():
            total = ensg_expr[gene][sample]
            for iso, val in iso_dict.items():
                rel = val / total if total > 0 else 0.0
                rel_expr[gene][iso][sample] = rel
                all_rel.append(rel)
                prev = rel_str[gene].get(iso, "")
                rel_str[gene][iso] = (prev + " " + f"{rel:.3f}") if prev else f"{rel:.3f}"

    if args.verbose:
        print(f"Processed expression: {len(ensg_expr)} genes, {len(all_rel)} isoform rel-values", file=sys.stderr)
    return ensg_expr, enst_expr, rel_expr, rel_str, all_rel, samples


def enrichment(data, min_exp, out_file=None):
    fh_out = gzip.open(out_file, 'wt') if out_file else None
    enrich = {}
    all_e   = []
    for gene, samp_dict in data.items():
        for sample, iso_dict in samp_dict.items():
            sorted_iso = sorted(iso_dict.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_iso) < 2:
                continue
            top_iso, top_val = sorted_iso[0]
            if top_val < min_exp:
                continue
            second_iso, sec_val = sorted_iso[1]
            enr = top_val / sec_val if sec_val > 0 else float('inf')
            if fh_out:
                fh_out.write(
                    f"{sample}\t{gene}\t{top_iso}\t{second_iso}\t"
                    f"{top_val:.3f}\t{sec_val:.3f}\t{enr:.3f}\n"
                )
            if enr >= args.minEnr:
                enrich.setdefault(gene, {}).setdefault(top_iso, {})[sample] = enr
                all_e.append(enr)
    if fh_out:
        fh_out.close()
    if args.verbose:
        print(f"Computed enrichment: {len(all_e)} events across {len(enrich)} genes", file=sys.stderr)
    return enrich, all_e


def main():
    sequences            = read_sequence_file(args.seq)
    miss_int, int_counts = read_isoform_int_file(args.iso, args.minString)
    canon_map            = read_canon_file(args.canon)

    ensg_pc, enst_pc, rel_pc, rel_str_pc, all_rel_pc, samples_pc = \
        read_expression_file(args.pcawg, sequences)
    ensg_gt, enst_gt, rel_gt, rel_str_gt, all_rel_gt, samples_gt = \
        read_expression_file(args.gtex, sequences)

    enrich_pc, _ = enrichment(enst_pc, args.minExp, args.pcawgEnrichmentFile)
    enrich_gt, _ = enrichment(enst_gt, args.minExp, args.gtexEnrichmentFile)

    output1 = []
    output2 = []
    pvals   = []

    for gene in sorted(enrich_pc):
        if gene not in rel_gt:
            continue
        for iso in sorted(enrich_pc[gene]):
            if iso not in rel_gt[gene]:
                continue
            if gene in enrich_gt and iso in enrich_gt[gene]:
                continue
            if gene not in enrich_gt or not enrich_gt[gene]:
                continue
            for sample, enr_pc in enrich_pc[gene][iso].items():
                vals_gt = list(rel_gt[gene][iso].values())
                if not vals_gt:
                    continue
                med_gt = statistics.median(vals_gt)
                rel_c  = rel_pc[gene][iso][sample]
                pos    = sum(1 for v in vals_gt if rel_c > v)
                neg    = sum(1 for v in vals_gt if rel_c < v)
                if pos + neg == 0:
                    continue
                pv = binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue

                # GTEx MDT analysis
                gtex_counts = {}
                mdi_vals    = []
                for niso in enrich_gt[gene]:
                    for ns in enrich_gt[gene][niso]:
                        n_enr = enrich_gt[gene][niso][ns]
                        if n_enr >= args.minEnr and sequences.get(niso) and sequences[niso] != sequences[iso]:
                            gtex_counts[niso] = gtex_counts.get(niso, 0) + 1
                            mdi_vals.append(n_enr)
                total_gt = len(vals_gt)
                skip     = True
                gtex_list= []
                for niso, cnt in sorted(gtex_counts.items(), key=lambda x: x[1], reverse=True):
                    if cnt >= total_gt * 0.5:
                        skip = False
                    gtex_list.append(f"{niso}:{cnt}")
                if skip:
                    continue
                gtexMDTs   = ",".join(gtex_list) if gtex_list else "-"
                med_mdi_gt = statistics.median(mdi_vals) if mdi_vals else "-"

                # Interaction disruption
                string_n     = int_counts.get(iso, 0)
                missed_can   = 0
                missed_com   = 0
                missed_list  = []
                for c in miss_int.get(iso, {}):
                    parts = c.split(':')
                    if len(parts) < 7:
                        continue
                    ensp2 = parts[5]
                    g2    = canon_map.get(ensp2)
                    if g2 in ensg_pc and sample in ensg_pc[g2]:
                        found = any(c in miss_int.get(niso, {}) for niso in enrich_gt[gene])
                        if found:
                            missed_com += 1
                        else:
                            missed_can += 1
                            missed_list.append(c)
                missed_str = ",".join(missed_list) if missed_list else "-"

                can_ensp = canon_map.get(iso, "-")
                num_gt   = len(vals_gt)

                o1 = (
                    f"{gene}\t{can_ensp}\t{iso}\t{sample}\t"
                    f"{gtexMDTs}\t{num_gt}\t{med_mdi_gt}\n"
                )
                o2 = (
                    f"{len(mdi_vals)}\t{enr_pc:.3f}\t{med_gt:.3f}\t"
                    f"{rel_c:.3f}\t{missed_str}\t{string_n}\t"
                    f"{missed_can}\t{missed_com}\n"
                )
                output1.append(o1)
                output2.append(o2)
                pvals.append(pv)

    if not pvals:
        if args.verbose:
            print("No P-values to adjust; exiting.", file=sys.stderr)
        sys.exit(0)

    reject, qvals, _, _ = multipletests(pvals, alpha=args.maxQ, method='fdr_bh')

    # Write output to gzipped file, one complete record per line
    out_file = 'interactionDisruptionInDominantTranscripts_min2.tsv.gz'
    with gzip.open(out_file, 'wt') as out_fh:
        header = [
            "#ENSG", "canENSP", "ENST", "Sample", "GTExMDTs", "N_GTEx", "MedianMDTenrGTEx",
            "Pvalue", "Qvalue", "N_GTeX_MDTs", "EnrPcawg", "MedRelExprGTEx",
            "RelExprPcawg", "MissedInts", "StringIntN", "MissedCanN", "MissedComN"
        ]
        out_fh.write("\t".join(header) + "\n")
        for o1_val, pv_val, qv_val, o2_val in zip(output1, pvals, qvals, output2):
            line = o1_val.rstrip("\n") + f"\t{pv_val:.3e}\t{qv_val:.3e}\t" + o2_val
            out_fh.write(line)

    if args.verbose:
        print(f"Output written to {out_file}", file=sys.stderr)


if __name__ == '__main__':
    main()