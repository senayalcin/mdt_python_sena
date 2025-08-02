import argparse
import gzip
from collections import defaultdict
import sys
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from statistics import median
import random


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Bu script PCAWG ve GTEx verilerini karşılaştırarak dominant transkriptleri analiz eder.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-pcawg", required=True, help="PCAWG RNA-seq TPM dosyası (örnek: pcawg.rnaseq.transcript.expr.tpm.tsv)")
    parser.add_argument("-gtex", required=True, help="GTEx RNA-seq TPM dosyası")
    parser.add_argument("-ensgFile", required=True, help="ENSG-ENSP-ENST eşleme dosyası")
    parser.add_argument("-canonFile", required=True, help="Canonical izoform bilgisi")
    parser.add_argument("-isoformIntFile", required=True, help="İzofrom etkileşim dosyası")
    parser.add_argument("-sequenceFile", required=True, help="ENSP sekans dosyası")
    parser.add_argument("-minStringScore", type=int, required=True, help="Minimum STRING skoru (örn: 900)")
    parser.add_argument("-enstColumnNuo", type=int, required=True, help="TPM dosyasında ENST kolon indexi (örn: 0)")
    parser.add_argument("-minEnrichment", type=float, required=True, help="Minimum zenginleşme katsayısı (örn: 2.0)")
    parser.add_argument("-maxQvalue", type=float, required=True, help="Maksimum q-değeri (örn: 0.01)")
    parser.add_argument("-minExpress", type=float, required=True, help="İfade eşik değeri (örn: 2.0)")

    parser.add_argument("-pvalueFile", help="İsteğe bağlı p-değeri dosyası")
    parser.add_argument("-randomPvalueFile", help="Rastgele p-değeri dosyası")
    parser.add_argument("-pcawgEnrichmentFile", help="PCAWG zenginleşme çıktısı")
    parser.add_argument("-gtexEnrichmentFile", help="GTEx zenginleşme çıktısı")
    parser.add_argument("-redundantFile", help="Redundant ENST dosyası")
    parser.add_argument("-verbose", action="store_true", help="Detaylı çıktı modu")

    return parser.parse_args()



def read_enst2ensg_file(enst2ensg_filepath):
    """
    ENST → ENSG eşleşmesini içeren bir sözlük döner.
    Dosya gzip formatındadır ve tab ile ayrılmıştır.
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
    TPM dosyasını okur ve:
    - gen düzeyi ifade
    - transkript düzeyi ifade
    - göreli ifadeler
    döndürür.
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

        # Redundant düzeltmesi
        if enst in redundant_dict:
            if verbose:
                print(f"Replacing redundant cDNA {enst} with {redundant_dict[enst]}")
            enst = redundant_dict[enst]

        # Format kontrolü
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
    Dominant transkript zenginleşme oranını hesaplar.

    Parametreler:
    - expression_data: dict[ensg][sample][enst] = expression (TPM)
    - min_expression: float, minimum MDT ifadesi (örn: 2.0)
    - min_enrichment: float, MDT'nin diğerine göre kaç kat fazla olması gerektiği
    - output_file_path: str, varsa sonuç dosyasına yazılır (gzipli .tsv.gz)

    Döner:
    - enrichment_dict: dict[ensg][dominant_enst][sample] = enrichment değeri
    - all_enrichment_values: list[float], istatistiksel testler için tüm değerler
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

            # dominant ve ikinci transkript
            dominant_enst, dominant_expr = sorted_ensts[0]
            other_enst, other_expr = sorted_ensts[1]

            # En az belirli bir ifade seviyesinde mi?
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
    Canonical protein ID'leri okur. Üç tür eşleştirme yapar:
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
    Her transkript için:
    - Kaybolan etkileşimleri (stringScore ≥ min)
    - Toplam etkileşim sayısını (kayıp + mevcut)

    döndürür.

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

            # Kayıp etkileşimler
            miss_interactions = miss_int_field.split(",")
            for int_str in miss_interactions:
                if ":" in int_str:
                    fields = int_str.split(":")
                    if len(fields) >= 8:
                        string_score = float(fields[6])
                        if string_score >= min_string_score:
                            miss_int_for_isoform[enst1][int_str] = 1
                            interaction_count[enst1] += 1

            # Mevcut etkileşimler
            exist_interactions = exist_int_field.split(",")
            for int_str in exist_interactions:
                if ":" in int_str:
                    fields = int_str.split(":")
                    if len(fields) >= 8:
                        string_score = float(fields[6])
                        if string_score >= min_string_score:
                            interaction_count[enst1] += 1

    return miss_int_for_isoform, interaction_count

def read_redundant_file(redundant_file_path, verbose=False):
    """
    Redundant transkriptleri eşleştirir.
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
    ENST → protein/aminoasit dizisi eşlemesini döndürür.
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
    Bir sayı listesinin medyanını döndürür.
    """
    sorted_vals = sorted(values)
    n = len(sorted_vals)

    if n == 0:
        return 0  # veya None?

    if n % 2 == 1:
        return sorted_vals[n // 2]
    else:
        return (sorted_vals[(n // 2) - 1] + sorted_vals[n // 2]) / 2



def main():
    args = parse_arguments()

    # Gerekli parametrelerin kontrolü (argparse zaten yapıyor ama buraya ekliyoruz gerekirse)
    required_params = [
        args.pcawg, args.gtex, args.ensgFile, args.canonFile,
        args.isoformIntFile, args.sequenceFile, args.enstColumnNuo,
        args.minStringScore, args.minEnrichment, args.maxQvalue, args.minExpress
    ]
    if not all(required_params):
        print("Eksik parametre var. Lütfen -h ile yardım alın.", file=sys.stderr)
        sys.exit(1)

    # Redundant ve eşleşme dosyaları
    enst2ensg = read_enst2ensg_file(args.ensgFile)
    redundant_map = read_redundant_file(args.redundantFile, verbose=args.verbose)
    sequences = read_sequence_file(args.sequenceFile, verbose=args.verbose)
    canonical = read_canon_file(args.canonFile)

    # İfade dosyaları
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

    # Etkileşim dosyası
    miss_iso, int_counts = read_isoform_interaction_file(args.isoformIntFile, args.minStringScore)

    # Enrichment hesaplama
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


###### burdan sonra toplu verdim kodu olmadı, sonra ayrı ayrı verip biz birleştirdik belki giritilerde sorun olabilir mantık hatası falanda ????
    output_part1 = []
    output_part2 = []
    pvalues_list = []

    for ensg in sorted(enrichment_pcawg):
        if ensg not in relstr_gtex:
            print(f"WARNING: {ensg} does not exist in GTEx data. Skipping.")
            continue

        for enst in sorted(enrichment_pcawg[ensg]):
            if enst not in relstr_gtex[ensg]:
                print(f"WARNING: {ensg} - {enst} not in GTEx data. Skipping.")
                continue

            for sample in sorted(enrichment_pcawg[ensg][enst]):
                enrichment_pcawg_val = enrichment_pcawg[ensg][enst][sample]
                rel_expr_pcawg = rel_pcawg[ensg][enst][sample]
                gtex_expr_str = relstr_gtex[ensg][enst]

                if enst in enrichment_gtex.get(ensg, {}):
                    continue  # GTEx’te MDT olarak geçiyorsa atla

                if not enrichment_gtex.get(ensg):
                    continue  # GTEx’te hiç MDT yoksa atla

                gtex_expr_values = list(map(float, gtex_expr_str.strip().split()))
                if not gtex_expr_values:
                    continue

                rel_expr_gtex_median = median(gtex_expr_values)
                pos = sum(1 for val in gtex_expr_values if rel_expr_pcawg > val)
                neg = sum(1 for val in gtex_expr_values if rel_expr_pcawg < val)

                # Binom testi
                pval = binomtest([pos, neg], p=0.5, alternative="two-sided").pvalue

                # Pvalue dosyasına yaz
                if pvalue_file:
                    pvalue_file.write(
                        f"{ensg}\t{enst}\t{sample}\t"
                        f"{rel_expr_gtex_median:.3f}\t{rel_expr_pcawg:.3f}\t"
                        f"{pos}\t{neg}\t{pval:.3e}\t{gtex_expr_str}\n"
                    )

                # Rastgele karşılaştırmalar
                for _ in range(100):
                    rand_enrichment = random.choice(all_enrich_pcawg)
                    rpos = sum(1 for val in gtex_expr_values if rand_enrichment > val)
                    rneg = sum(1 for val in gtex_expr_values if rand_enrichment < val)
                    rand_pval = binomtest([rpos, rneg], p=0.5).pvalue
                    if random_pvalue_file:
                        random_pvalue_file.write(
                            f"{ensg}\t{enst}\t{sample}\t"
                            f"{rel_expr_gtex_median:.3f}\t{rand_enrichment:.3f}\t"
                            f"{rpos}\t{rneg}\t{rand_pval:.3e}\n"
                        )

                if pval >= args.maxQvalue:
                    continue

                # GTEx MDT'lerini bul
                gtex_mdts = {}
                mdi_enrichments = []
                if ensg in enrichment_gtex:
                    for n_enst in enrichment_gtex[ensg]:
                        for n_sample, enrich_val in enrichment_gtex[ensg][n_enst].items():
                            if (
                                    enrich_val >= args.minEnrichment
                                    and sequences.get(n_enst) != sequences.get(enst)
                            ):
                                gtex_mdts[n_enst] = gtex_mdts.get(n_enst, 0) + 1
                                mdi_enrichments.append(enrich_val)

                # %50 kuralı
                total_normals = len(gtex_expr_values)
                skip = True
                gtex_mdts_str = ""
                for g in sorted(gtex_mdts, key=gtex_mdts.get, reverse=True):
                    gtex_mdts_str += f"{g}:{gtex_mdts[g]},"
                    if gtex_mdts[g] >= 0.5 * total_normals:
                        skip = False
                if skip:
                    continue
                gtex_mdts_str = gtex_mdts_str.rstrip(",")

                # MDT’ler için median enrichment
                median_gtex_enrichment = "-"
                if mdi_enrichments:
                    median_gtex_enrichment = f"{median(mdi_enrichments):.3f}"

                # interaction disruption kontrol
                missed_ints = []
                missed_common = 0
                missed_cancer = 0
                string_int_n = int_counts.get(enst, -1)

                if enst in miss_iso:
                    for cint in miss_iso[enst]:
                        parts = cint.split(":")
                        if len(parts) < 8:
                            continue
                        _, _, _, _, _, ensp2, _, _ = parts
                        partner_ensg = canonical.get(ensp2)
                        if partner_ensg and sample in ensg_pcawg.get(partner_ensg, {}):
                            found = False
                            if ensg in enrichment_gtex:
                                for alt_enst in enrichment_gtex[ensg]:
                                    for alt_sample in enrichment_gtex[ensg][alt_enst]:
                                        if (
                                                enrichment_gtex[ensg][alt_enst][alt_sample] >= args.minEnrichment
                                                and cint in miss_iso.get(alt_enst, {})
                                        ):
                                            found = True
                            if found:
                                missed_common += 1
                            else:
                                missed_ints.append(cint)
                                missed_cancer += 1

                missed_str = ",".join(missed_ints) if missed_ints else "-"

                # Veri dizilerine ekle
                output_part1.append(
                    f"{ensg}\t{canonical.get(enst, '-')}\t{enst}\t{sample}\t"
                    f"{gtex_mdts_str or '-'}\t{total_normals}\t{median_gtex_enrichment}"
                )
                output_part2.append(
                    f"{len(mdi_enrichments)}\t{enrichment_pcawg_val:.3f}\t"
                    f"{rel_expr_gtex_median:.3f}\t{rel_expr_pcawg:.3f}\t"
                    f"{missed_str}\t{string_int_n}\t{missed_cancer}\t{missed_common}"
                )
                pvalues_list.append(pval)

    # Q-value düzeltme ve çıktı
    if pvalues_list:
        _, qvalues, _, _ = multipletests(pvalues_list, alpha=args.maxQvalue, method='fdr_bh')
        for i in range(len(output_part1)):
            if qvalues[i] < args.maxQvalue:
                print(f"{output_part1[i]}\t{pvalues_list[i]:.3e}\t{qvalues[i]:.3e}\t{output_part2[i]}")

    if pvalue_file:
        pvalue_file.close()
    if random_pvalue_file:
        random_pvalue_file.close()







