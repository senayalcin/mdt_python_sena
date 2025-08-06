#!/usr/bin/env python3
"""
Author: Abdullah Kahraman (converted to Python)
Date: 06.02.2018

Determines cancer-specific MDT based on 2-fold expression difference to
minor dominant transcripts and subsequently compares the relative
expression difference to expression differences in GTEx using Sign-test.
"""

import sys
import os
import argparse
import gzip
import random
import statistics
from scipy import stats
from scipy.stats import binomtest  # R binom.test ile uyumlu
import numpy as np
from statsmodels.stats.multitest import fdrcorrection


class MDTAnalyzer:
    def __init__(self):
        self.verbose = False
        self.enst_column_no = 1
        self.min_string_score = 900
        self.min_enrichment = 2.0
        self.max_qvalue = 0.01
        self.min_express = 2.0

        # File paths
        self.enst2ensg_file = None
        self.express_pcawg_file = None
        self.express_gtex_file = None
        self.isoform_int_file = None
        self.canon_file = None
        self.sequence_file = None
        self.pvalue_file = None
        self.random_pvalue_file = None
        self.pcawg_enrichment_file = None
        self.gtex_enrichment_file = None
        self.redundant_file = None

    def parse_arguments(self):
        parser = argparse.ArgumentParser(description='Most Dominant Transcript Switch Analysis')

        parser.add_argument('--pcawg', required=True,
                            help='PCAWG expression file (e.g. pcawg.rnaseq.transcript.expr.tpm.tsv.gz)')
        parser.add_argument('--gtex', required=True,
                            help='GTEx expression file (e.g. GTEX_v4.pcawg.transcripts.tpm.tsv.gz)')
        parser.add_argument('--ensgFile', required=True,
                            help='Ensembl ID mapping file (e.g. ensg_ensp_enst_ense_geneName_v75.tsv.gz)')
        parser.add_argument('--canonFile', required=True,
                            help='Canonical isoforms file (e.g. canonical_isoforms_ensg_ensp_enst_geneName_v75.tsv.gz)')
        parser.add_argument('--isoformIntFile', required=True,
                            help='Isoform interaction file (e.g. interactionsInIsoforms_900_2.tsv.gz)')
        parser.add_argument('--sequenceFile', required=True,
                            help='Sequence file (e.g. ensp_ensg_enst_sequence.tsv.gz)')
        parser.add_argument('--minStringScore', type=int, required=True,
                            help='Minimum STRING score (e.g. 900)')
        parser.add_argument('--enstColumnNo', type=int, default=1,
                            help='Column number of ENST in expression files (starts with 0, default=1)')
        parser.add_argument('--minEnrichment', type=float, required=True,
                            help='Minimum fold difference (e.g. 2)')
        parser.add_argument('--maxQvalue', type=float, required=True,
                            help='Maximum Q-value/P-value (e.g. 0.05)')
        parser.add_argument('--minExpress', type=float, required=True,
                            help='Minimum expression value (e.g. 1)')
        parser.add_argument('--pvalueFile',
                            help='P-values output file (e.g. pvalues.tsv)')
        parser.add_argument('--randomPvalueFile',
                            help='Random p-values output file (e.g. random_pvalues.tsv)')
        parser.add_argument('--pcawgEnrichmentFile',
                            help='PCAWG enrichment output file (e.g. pcawg_express.tsv.gz)')
        parser.add_argument('--gtexEnrichmentFile',
                            help='GTEx enrichment output file (e.g. gtex_express.tsv.gz)')
        parser.add_argument('--redundantFile',
                            help='Redundant ENST file (e.g. redundantENST.txt)')
        parser.add_argument('--verbose', action='store_true',
                            help='Print additional information')

        args = parser.parse_args()

        # Set attributes from arguments
        self.express_pcawg_file = args.pcawg
        self.express_gtex_file = args.gtex
        self.enst2ensg_file = args.ensgFile
        self.canon_file = args.canonFile
        self.isoform_int_file = args.isoformIntFile
        self.sequence_file = args.sequenceFile
        self.min_string_score = args.minStringScore
        self.enst_column_no = args.enstColumnNo
        self.min_enrichment = args.minEnrichment
        self.max_qvalue = args.maxQvalue
        self.min_express = args.minExpress
        self.pvalue_file = args.pvalueFile
        self.random_pvalue_file = args.randomPvalueFile
        self.pcawg_enrichment_file = args.pcawgEnrichmentFile
        self.gtex_enrichment_file = args.gtexEnrichmentFile
        self.redundant_file = args.redundantFile
        self.verbose = args.verbose

    def read_enst2ensg_file(self):
        """Read ENST to ENSG mapping file"""
        enst2ensg = {}

        if self.enst2ensg_file.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(self.enst2ensg_file, 'rt') as f:
            for line in f:
                if not line.startswith('ENS'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    enst2ensg[parts[2]] = parts[0]

        return enst2ensg

    def read_redundant_file(self):
        """Read redundant ENST file"""
        redundant = {}
        if not self.redundant_file:
            return redundant

        with open(self.redundant_file, 'r') as f:
            for line in f:
                if not line.startswith('EN'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    redundant[parts[2]] = parts[1]

        if self.verbose:
            print(f"Found {len(redundant)} redundant ENSTs", file=sys.stderr)
        return redundant

    def read_expression_file(self, filename):
        """Read expression file and calculate relative expressions"""
        enst2ensg = self.read_enst2ensg_file()
        redundant = self.read_redundant_file()

        enst_express = {}
        ensg_express = {}

        if filename.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(filename, 'rt') as f:
            samples = []
            n = 0

            for line in f:
                line = line.strip()

                # Read header with sample names
                if n == 0:
                    samples = line.split('\t')
                    n += 1
                    continue

                parts = line.split('\t')
                enst = parts[self.enst_column_no]

                # Handle redundant transcripts
                if enst in redundant:
                    if self.verbose:
                        print(f"Replacing redundant cDNA {enst} with {redundant[enst]}")
                    enst = redundant[enst]

                if not enst.startswith('ENST0'):
                    continue

                # Remove version number
                enst = enst.split('.')[0]

                if enst not in enst2ensg:
                    print(f"WARNING: {enst} does not exist in {self.enst2ensg_file}. Skipping transcript.",
                          file=sys.stderr)
                    continue

                ensg = enst2ensg[enst]

                # Process expression values starting from column after ENST
                for i in range(self.enst_column_no + 1, len(parts)):
                    if i >= len(samples):
                        break

                    expression = parts[i]
                    if expression == "NA":
                        expression = 0.0
                    else:
                        try:
                            expression = float(expression)
                        except ValueError:
                            expression = 0.0

                    sample = samples[i]

                    # Initialize nested dictionaries
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

        # Calculate relative transcript expression - PERL GİBİ SIRALI
        rel_enst_express = {}
        rel_enst_express_string = {}
        all_rel_enst_express = []

        # KRITIK: Perl scriptindeki gibi sıralama yapılmalı
        for ensg in sorted(enst_express.keys()):
            for sample in sorted(enst_express[ensg].keys()):
                gene_express = ensg_express[ensg][sample]

                for enst in sorted(enst_express[ensg][sample].keys()):
                    trans_express = enst_express[ensg][sample][enst]

                    rel_express = 0.0
                    if gene_express > 0:
                        rel_express = trans_express / gene_express

                    # Initialize nested dictionaries
                    if ensg not in rel_enst_express:
                        rel_enst_express[ensg] = {}
                    if enst not in rel_enst_express[ensg]:
                        rel_enst_express[ensg][enst] = {}

                    rel_enst_express[ensg][enst][sample] = rel_express
                    all_rel_enst_express.append(rel_express)

                    # Build string representation - PERL GİBİ
                    if ensg not in rel_enst_express_string:
                        rel_enst_express_string[ensg] = {}

                    if enst in rel_enst_express_string[ensg]:
                        rel_enst_express_string[ensg][enst] += f" {rel_express:.3f}"
                    else:
                        rel_enst_express_string[ensg][enst] = f"{rel_express:.3f}"

        return (ensg_express, enst_express, rel_enst_express,
                rel_enst_express_string, all_rel_enst_express)

    def calculate_enrichment(self, express, minimum_express, output_file=None):
        """Calculate enrichment values"""
        if output_file:
            if output_file.endswith('.gz'):
                out_file = gzip.open(output_file, 'wt')
            else:
                out_file = open(output_file, 'w')

        enrichment = {}
        all_enrichment = []

        for ensg in sorted(express.keys()):
            for sample in sorted(express[ensg].keys()):
                most_dominant_enst = "-"
                most_dominant_enst_express = 0
                enrich = -1
                no_enst = len(express[ensg][sample])

                # Sort isoforms by highest expression
                sorted_transcripts = sorted(express[ensg][sample].items(),
                                            key=lambda x: x[1], reverse=True)

                n = 0
                for enst, enst_express in sorted_transcripts:
                    n += 1

                    # First transcript is the one with highest expression (potential MDT)
                    if n == 1:
                        most_dominant_enst = enst
                        most_dominant_enst_express = enst_express
                        continue

                    enrichment_val = most_dominant_enst_express

                    # Ignore transcripts with insignificant expression
                    if most_dominant_enst_express >= minimum_express:
                        if enst_express > 0:
                            enrichment_val = most_dominant_enst_express / enst_express

                        if output_file:
                            out_file.write(f"{sample}\t{ensg}\t{most_dominant_enst}\t{enst}\t"
                                           f"{most_dominant_enst_express}\t{enst_express}\t{enrichment_val}\n")

                        if enrichment_val >= self.min_enrichment:
                            if ensg not in enrichment:
                                enrichment[ensg] = {}
                            if most_dominant_enst not in enrichment[ensg]:
                                enrichment[ensg][most_dominant_enst] = {}

                            enrichment[ensg][most_dominant_enst][sample] = enrichment_val
                            all_enrichment.append(enrichment_val)

                    if n == 2:  # FIXED: Only compare with second highest
                        break

        if output_file:
            out_file.close()

        return enrichment, all_enrichment

    def read_canon_file(self):
        """Read canonical isoforms file"""
        canon = {}

        if self.canon_file.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(self.canon_file, 'rt') as f:
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

    def read_isoform_int_file(self):
        """Read isoform interaction file"""
        miss_int_4_isoform = {}
        int_n = {}

        if self.isoform_int_file.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(self.isoform_int_file, 'rt') as f:
            for line in f:
                if not line.startswith('EN'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 10:
                    (ensg1, ensp1, enst1, string_ensp1, string_int_n,
                     exist_int_n, miss_int_n, rel_miss_int_n, exist_int, miss_int) = parts[:10]

                    int_n[enst1] = 0

                    # Process missed interactions
                    if miss_int:
                        for interaction in miss_int.split(','):
                            if ':' in interaction:
                                int_parts = interaction.split(':')
                                if len(int_parts) >= 8:
                                    string_score = int(int_parts[6])
                                    if string_score >= self.min_string_score:
                                        if enst1 not in miss_int_4_isoform:
                                            miss_int_4_isoform[enst1] = {}
                                        miss_int_4_isoform[enst1][interaction] = 1
                                        int_n[enst1] += 1

                    # Process existing interactions
                    if exist_int:
                        for interaction in exist_int.split(','):
                            if ':' in interaction:
                                int_parts = interaction.split(':')
                                if len(int_parts) >= 8:
                                    string_score = int(int_parts[6])
                                    if string_score >= self.min_string_score:
                                        int_n[enst1] += 1

        return miss_int_4_isoform, int_n

    def read_sequence_file(self):
        """Read sequence file"""
        sequences = {}

        if self.sequence_file.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        with opener(self.sequence_file, 'rt') as f:
            for line in f:
                if not line.startswith('EN'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    ensp, ensg, enst, sequence = parts[:4]
                    sequences[enst] = sequence

        if self.verbose:
            print(f"Found {len(sequences)} sequences", file=sys.stderr)

        return sequences

    def median(self, values):
        """Calculate median of values"""
        if not values:
            return 0
        return statistics.median(values)

    def run_analysis(self):
        """Main analysis function"""
        # Set random seed
        random.seed(23)

        if self.verbose:
            print(f"Reading {self.sequence_file} ... ", end="", file=sys.stderr)
        sequences = self.read_sequence_file()
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print(f"Reading in {self.express_pcawg_file} ... ", end="", file=sys.stderr)
        (ensg_express_pcawg, enst_express_pcawg, rel_enst_express_pcawg,
         rel_enst_express_string_pcawg, all_rel_enst_express_pcawg) = self.read_expression_file(self.express_pcawg_file)
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print(f"Reading in {self.express_gtex_file} ... ", end="", file=sys.stderr)
        (ensg_express_gtex, enst_express_gtex, rel_enst_express_gtex,
         rel_enst_express_string_gtex, all_rel_enst_express_gtex) = self.read_expression_file(self.express_gtex_file)
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print(f"Reading {self.canon_file} ... ", end="", file=sys.stderr)
        canonical = self.read_canon_file()
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print(f"Reading {self.isoform_int_file} ... ", end="", file=sys.stderr)
        miss_iso_int, int_n = self.read_isoform_int_file()
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print("Calculating Enrichment for PCAWG ... ", end="", file=sys.stderr)
        enrichment_pcawg, all_enrichment_pcawg = self.calculate_enrichment(
            enst_express_pcawg, self.min_express, self.pcawg_enrichment_file)
        if self.verbose:
            print("done", file=sys.stderr)

        if self.verbose:
            print("Calculating Enrichment for GTEx ... ", end="", file=sys.stderr)
        enrichment_gtex, all_enrichment_gtex = self.calculate_enrichment(
            enst_express_gtex, self.min_express * 0.1, self.gtex_enrichment_file)
        if self.verbose:
            print("done", file=sys.stderr)

        # Print file header
        print("#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\t"
              "GTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\tNumberOfGTExMDT\tcMDTenrichment\t"
              "GTExMDTmedianRelExpression\tcMDTrelExpression\t"
              "cMDTuniqMissedInt\tTotalNumberOfSTRINGint\t"
              "cMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt")

        # Open output files if specified
        pvalue_file_handle = None
        random_pvalue_file_handle = None

        if self.pvalue_file:
            pvalue_file_handle = open(self.pvalue_file, 'w')
            pvalue_file_handle.write(
                "#ENSG\tcMDT\ttRNAseqAliquotID\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
                "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\tRelMDTexpressionInGTEx\n")

        if self.random_pvalue_file:
            random_pvalue_file_handle = open(self.random_pvalue_file, 'w')
            random_pvalue_file_handle.write(
                "#ENSG\tDomCancerTrans\tCancerSampleId\tMedianRelMDTexpressionInNormal\tRelMDTexpressionInCancer\t"
                "PCAWGhigherExpressionN\tPCAWGlowerExpressionN\tPvalue\n")

        output_part1 = []
        output_part2 = []
        pvalues_list = []

        # Main analysis loop - PERL GİBİ SIRALAMA İÇİN
        for ensg in sorted(enrichment_pcawg.keys()):
            if ensg not in rel_enst_express_string_gtex:
                print(f"WARNING: {ensg} does not exists in {self.express_gtex_file}. Skipping {ensg}.",
                      file=sys.stderr)
                continue

            for enst in sorted(enrichment_pcawg[ensg].keys()):
                if enst not in rel_enst_express_string_gtex[ensg]:
                    print(f"WARNING: {ensg} - {enst} does not exists in {self.express_gtex_file}. Skipping it.",
                          file=sys.stderr)
                    continue

                # KRITIK: Sample'ları da sıralı olarak işle
                for sample in sorted(enrichment_pcawg[ensg][enst].keys()):
                    enrichment_pcawg_val = enrichment_pcawg[ensg][enst][sample]
                    enst_rel_express_pcawg = rel_enst_express_pcawg[ensg][enst][sample]
                    enst_rel_expresses_gtex = rel_enst_express_string_gtex[ensg][enst]

                    # FIXED: Check conditions like in Perl
                    if (enst not in enrichment_gtex.get(ensg, {}) and
                            ensg in enrichment_gtex and len(enrichment_gtex[ensg]) > 0):

                        enst_rel_expresses_gtex_list = [float(x) for x in enst_rel_expresses_gtex.split()]

                        rel_expression_gtex_median = 0
                        if enst_rel_expresses_gtex_list:
                            rel_expression_gtex_median = self.median(enst_rel_expresses_gtex_list)

                        pos = sum(1 for x in enst_rel_expresses_gtex_list if enst_rel_express_pcawg > x)
                        neg = sum(1 for x in enst_rel_expresses_gtex_list if enst_rel_express_pcawg < x)

                        # Binomial test - PERL R İLE AYNI SONUÇ İÇİN
                        if pos + neg > 0:
                            # R'daki binom.test ile aynı sonuç için scipy.stats.binom_test kullan
                            from scipy.stats import binomtest
                            try:
                                pvalue = binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue
                            except:
                                # Fallback to binomtest if binom_test not available
                                pvalue = stats.binomtest(pos, pos + neg, p=0.5, alternative='two-sided').pvalue
                        else:
                            pvalue = 1.0

                        if pvalue_file_handle:
                            pvalue_file_handle.write(f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t"
                                                     f"{enst_rel_express_pcawg:.3f}\t{pos}\t{neg}\t{pvalue:.3e}\t{enst_rel_expresses_gtex}\n")

                        # Random permutation tests
                        if random_pvalue_file_handle:
                            for j in range(100):
                                random_enrichment_pcawg = random.choice(all_enrichment_pcawg)

                                r_pos = sum(1 for x in enst_rel_expresses_gtex_list if random_enrichment_pcawg > x)
                                r_neg = sum(1 for x in enst_rel_expresses_gtex_list if random_enrichment_pcawg < x)

                                if r_pos + r_neg > 0:
                                    from scipy.stats import binomtest
                                    try:
                                        random_pvalue = binomtest(r_pos, r_pos + r_neg, p=0.5, alternative='two-sided').pvalue
                                    except:
                                        random_pvalue = stats.binomtest(r_pos, r_pos + r_neg, p=0.5, alternative='two-sided').pvalue
                                else:
                                    random_pvalue = 1.0

                                random_pvalue_file_handle.write(
                                    f"{ensg}\t{enst}\t{sample}\t{rel_expression_gtex_median:.3f}\t"
                                    f"{random_enrichment_pcawg:.3f}\t{r_pos}\t{r_neg}\t{random_pvalue:.3e}\n")

                        # Continue with significant p-values
                        if pvalue < self.max_qvalue:
                            pvalues_list.append(pvalue)

                            # Process missed interactions
                            missed_int = ""
                            string_int_n = -1
                            if enst in int_n:
                                string_int_n = int_n[enst]

                            missed_cancer_int_n = -1
                            missed_common_int_n = -1
                            gtex_mdts = {}
                            mdi_enrichments_gtex = []

                            # Check GTEx MDTs
                            if ensg in enrichment_gtex:
                                for n_enst in enrichment_gtex[ensg]:
                                    for n_sample in enrichment_gtex[ensg][n_enst]:
                                        n_enrichment_gtex = enrichment_gtex[ensg][n_enst][n_sample]
                                        if (n_enrichment_gtex >= self.min_enrichment and
                                                n_enst in sequences and enst in sequences and
                                                sequences[n_enst] != sequences[enst]):
                                            gtex_mdts[n_enst] = gtex_mdts.get(n_enst, 0) + 1
                                            mdi_enrichments_gtex.append(n_enrichment_gtex)

                            gtex_mdts_str = ""
                            skip = True
                            for g in sorted(gtex_mdts.keys(), key=lambda x: gtex_mdts[x], reverse=True):
                                gtex_mdts_str += f"{g}:{gtex_mdts[g]},"
                                if gtex_mdts[g] >= len(enst_rel_expresses_gtex_list) * 0.5:
                                    skip = False

                            if skip:
                                continue

                            gtex_mdts_str = gtex_mdts_str.rstrip(',')
                            if gtex_mdts_str == "":
                                gtex_mdts_str = "-"

                            median_mdt_enrichment_gtex = "-"
                            if mdi_enrichments_gtex:
                                median_mdt_enrichment_gtex = f"{self.median(mdi_enrichments_gtex):.3f}"

                            # Check interactions
                            if enst in miss_iso_int:
                                missed_cancer_int_n = 0
                                missed_common_int_n = 0

                                for c_int in miss_iso_int[enst]:
                                    int_parts = c_int.split(':')
                                    if len(int_parts) >= 8:
                                        ensp2 = int_parts[5]

                                        if (ensp2 in canonical and canonical[ensp2] in ensg_express_pcawg and
                                                sample in ensg_express_pcawg[canonical[ensp2]]):

                                            found = False
                                            if ensg in enrichment_gtex:
                                                for n_enst in enrichment_gtex[ensg]:
                                                    for n_sample in enrichment_gtex[ensg][n_enst]:
                                                        enrichment_gtex_val = enrichment_gtex[ensg][n_enst][n_sample]
                                                        if (enrichment_gtex_val >= self.min_enrichment and
                                                                n_enst in miss_iso_int and c_int in miss_iso_int[
                                                                    n_enst]):
                                                            found = True
                                                            break
                                                    if found:
                                                        break

                                            if found:
                                                missed_common_int_n += 1
                                            else:
                                                missed_int += f"{c_int},"
                                                missed_cancer_int_n += 1

                            missed_int = missed_int.rstrip(',')
                            if missed_int == "":
                                missed_int = "-"

                            # Store output parts
                            canon_enst = canonical.get(enst, "-")
                            output_part1.append(f"{ensg}\t{canon_enst}\t{enst}\t{sample}\t"
                                                f"{gtex_mdts_str}\t{len(enst_rel_expresses_gtex_list)}\t{median_mdt_enrichment_gtex}")

                            output_part2.append(f"{len(mdi_enrichments_gtex)}\t{enrichment_pcawg_val:.3f}\t"
                                                f"{rel_expression_gtex_median:.3f}\t{enst_rel_express_pcawg:.3f}\t"
                                                f"{missed_int}\t{string_int_n}\t"
                                                f"{missed_cancer_int_n}\t{missed_common_int_n}")

        # FDR correction - FIXED: Use BH method like in Perl
        if pvalues_list:
            # Using Benjamini-Hochberg method like the Perl R p.adjust function
            rejected, qvalues_corrected = fdrcorrection(pvalues_list, alpha=1.0, method='indep')
        else:
            qvalues_corrected = []

        # Output results
        for i in range(len(output_part1)):
            if i < len(qvalues_corrected) and qvalues_corrected[i] < self.max_qvalue:
                print(f"{output_part1[i]}\t{pvalues_list[i]:.3e}\t{qvalues_corrected[i]:.3e}\t{output_part2[i]}")

        # Close output files
        if pvalue_file_handle:
            pvalue_file_handle.close()
        if random_pvalue_file_handle:
            random_pvalue_file_handle.close()


def main():
    """Main function"""
    analyzer = MDTAnalyzer()
    analyzer.parse_arguments()
    analyzer.run_analysis()


if __name__ == "__main__":
    main()