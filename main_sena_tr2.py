#!/usr/bin/env python3
"""
PCAWG ve GTEx RNA-seq verilerini karşılaştırarak dominant transkriptleri analiz eden script.

Bu script:
- PCAWG (kanser) ve GTEx (normal) RNA-seq TPM verilerini okur
- Dominant transkriptleri (MDT) belirler
- Zenginleşme analizleri yapar
- İstatistiksel testlerle anlamlı farkları bulur
- Protein etkileşim bozulmalarını analiz eder
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median
import random
import logging

# Bilimsel hesaplamalar için
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

# Logging konfigürasyonu
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_arguments():
    """Komut satırı argümanlarını parse et."""
    parser = argparse.ArgumentParser(
        description="PCAWG ve GTEx verilerini karşılaştırarak dominant transkriptleri analiz eder.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Zorunlu dosya parametreleri
    required_group = parser.add_argument_group('Zorunlu dosyalar')
    required_group.add_argument("-pcawg", required=True, type=Path,
                                help="PCAWG RNA-seq TPM dosyası (gzip)")
    required_group.add_argument("-gtex", required=True, type=Path,
                                help="GTEx RNA-seq TPM dosyası (gzip)")
    required_group.add_argument("-ensgFile", required=True, type=Path,
                                help="ENSG-ENSP-ENST eşleme dosyası (gzip)")
    required_group.add_argument("-canonFile", required=True, type=Path,
                                help="Canonical izoform bilgi dosyası (gzip)")
    required_group.add_argument("-isoformIntFile", required=True, type=Path,
                                help="İzoform etkileşim dosyası (gzip)")
    required_group.add_argument("-sequenceFile", required=True, type=Path,
                                help="ENSP protein sekans dosyası (gzip)")

    # Zorunlu analiz parametreleri
    analysis_group = parser.add_argument_group('Analiz parametreleri')
    analysis_group.add_argument("-minStringScore", type=int, required=True,
                                help="Minimum STRING skoru (örn: 900)")
    analysis_group.add_argument("-enstColumnNuo", type=int, required=True,
                                help="TPM dosyasında ENST kolon indexi (genellikle 0)")
    analysis_group.add_argument("-minEnrichment", type=float, required=True,
                                help="Minimum zenginleşme katsayısı (örn: 2.0)")
    analysis_group.add_argument("-maxQvalue", type=float, required=True,
                                help="Maksimum q-değeri eşiği (örn: 0.01)")
    analysis_group.add_argument("-minExpress", type=float, required=True,
                                help="Minimum ifade eşik değeri (TPM) (örn: 2.0)")

    # İsteğe bağlı çıktı dosyaları
    output_group = parser.add_argument_group('Çıktı dosyaları (isteğe bağlı)')
    output_group.add_argument("-pvalueFile", type=Path,
                              help="P-değerleri çıktı dosyası")
    output_group.add_argument("-randomPvalueFile", type=Path,
                              help="Rastgele karşılaştırma p-değerleri dosyası")
    output_group.add_argument("-pcawgEnrichmentFile", type=Path,
                              help="PCAWG zenginleşme detay çıktısı")
    output_group.add_argument("-gtexEnrichmentFile", type=Path,
                              help="GTEx zenginleşme detay çıktısı")
    output_group.add_argument("-redundantFile", type=Path,
                              help="Redundant ENST eşleme dosyası")

    # Diğer seçenekler
    parser.add_argument("-verbose", action="store_true",
                        help="Detaylı log çıktıları")

    return parser.parse_args()


def validate_files(args):
    """Dosya varlığını kontrol et."""
    required_files = [
        args.pcawg, args.gtex, args.ensgFile, args.canonFile,
        args.isoformIntFile, args.sequenceFile
    ]

    optional_files = [args.redundantFile] if args.redundantFile else []

    for file_path in required_files:
        if not file_path.exists():
            raise FileNotFoundError(f"Zorunlu dosya bulunamadı: {file_path}")

    for file_path in optional_files:
        if not file_path.exists():
            logger.warning(f"İsteğe bağlı dosya bulunamadı: {file_path}")


def read_enst2ensg_mapping(filepath, verbose=False):
    """ENST → ENSG eşleşmelerini oku."""
    enst2ensg = {}

    if verbose:
        logger.info(f"ENST-ENSG eşleme dosyası okunuyor: {filepath}")

    try:
        with gzip.open(filepath, "rt") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or not line.startswith("ENS"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 3:
                    ensg_id = parts[0]
                    enst_id = parts[2]
                    enst2ensg[enst_id] = ensg_id

        if verbose:
            logger.info(f"Toplam {len(enst2ensg)} ENST-ENSG eşleşmesi okundu")

        return enst2ensg

    except Exception as e:
        logger.error(f"ENST-ENSG dosyası okuma hatası: {e}")
        raise


def read_redundant_mapping(filepath, verbose=False):
    """Redundant transkript eşleşmelerini oku."""
    if not filepath:
        return {}

    redundant_map = {}

    if verbose:
        logger.info(f"Redundant eşleme dosyası okunuyor: {filepath}")

    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith("EN"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 3:
                    ensg, enst_selected, enst_redundant = parts[:3]
                    redundant_map[enst_redundant] = enst_selected

        if verbose:
            logger.info(f"Toplam {len(redundant_map)} redundant eşleşme okundu")

        return redundant_map

    except Exception as e:
        logger.error(f"Redundant dosyası okuma hatası: {e}")
        return {}


def read_expression_data(filepath, enst2ensg, redundant_map, enst_col_index=0, verbose=False):
    """TPM ifade verilerini oku ve işle."""
    if verbose:
        logger.info(f"İfade dosyası okunuyor: {filepath}")

    enst_expression = defaultdict(lambda: defaultdict(dict))
    ensg_expression = defaultdict(lambda: defaultdict(float))

    try:
        with gzip.open(filepath, "rt") as f:
            header = f.readline().strip().split("\t")
            sample_names = header[enst_col_index + 1:]

            for line_num, line in enumerate(f, 2):
                parts = line.strip().split("\t")
                if len(parts) <= enst_col_index:
                    continue

                enst = parts[enst_col_index]

                # Redundant düzeltmesi
                if enst in redundant_map:
                    if verbose and line_num <= 10:  # İlk birkaç örneği göster
                        logger.debug(f"Redundant ENST değiştiriliyor: {enst} -> {redundant_map[enst]}")
                    enst = redundant_map[enst]

                # Format kontrolü ve temizleme
                if "ENST0" not in enst:
                    continue
                enst = enst.split(".")[0]  # Version numarasını kaldır

                if enst not in enst2ensg:
                    if verbose:
                        logger.debug(f"ENST-ENSG eşleşmesi bulunamadı: {enst}")
                    continue

                ensg = enst2ensg[enst]

                # TPM değerlerini işle
                for i, tpm_str in enumerate(parts[enst_col_index + 1:]):
                    if i >= len(sample_names):
                        break

                    sample = sample_names[i]
                    tpm_value = 0.0 if tpm_str == "NA" else float(tpm_str)

                    enst_expression[ensg][sample][enst] = enst_expression[ensg][sample].get(enst, 0) + tpm_value
                    ensg_expression[ensg][sample] += tpm_value

        # Göreceli ifadeleri hesapla
        rel_expression = defaultdict(lambda: defaultdict(dict))
        rel_expression_strings = defaultdict(lambda: defaultdict(str))
        all_relative_values = []

        for ensg in enst_expression:
            for sample in enst_expression[ensg]:
                gene_total = ensg_expression[ensg][sample]

                for enst in enst_expression[ensg][sample]:
                    transcript_expr = enst_expression[ensg][sample][enst]
                    rel_value = transcript_expr / gene_total if gene_total > 0 else 0.0

                    rel_expression[ensg][enst][sample] = rel_value
                    all_relative_values.append(rel_value)

                    # String formatında da sakla
                    if rel_expression_strings[ensg][enst]:
                        rel_expression_strings[ensg][enst] += f" {rel_value:.3f}"
                    else:
                        rel_expression_strings[ensg][enst] = f"{rel_value:.3f}"

        if verbose:
            logger.info(
                f"Toplam {len(ensg_expression)} gen ve {sum(len(samples) for samples in enst_expression.values())} örnek işlendi")

        return ensg_expression, enst_expression, rel_expression, rel_expression_strings, all_relative_values

    except Exception as e:
        logger.error(f"İfade dosyası okuma hatası: {e}")
        raise


def compute_enrichment_analysis(expression_data, min_expression, min_enrichment, output_file=None, verbose=False):
    """Dominant transkript zenginleşme analizi yap."""
    enrichment_results = defaultdict(lambda: defaultdict(dict))
    all_enrichment_values = []

    if verbose:
        logger.info("Zenginleşme analizi başlatılıyor...")

    output_handle = None
    if output_file:
        output_handle = gzip.open(output_file, "wt")
        output_handle.write("Sample\tENSG\tDominantENST\tSecondENST\tDominantTPM\tSecondTPM\tEnrichment\n")

    try:
        processed_genes = 0
        for ensg in sorted(expression_data):
            for sample in sorted(expression_data[ensg]):
                transcripts = expression_data[ensg][sample]

                # En yüksek ifadeli transkriptleri sırala
                sorted_transcripts = sorted(transcripts.items(), key=lambda x: x[1], reverse=True)

                if len(sorted_transcripts) < 2:
                    continue

                dominant_enst, dominant_tpm = sorted_transcripts[0]
                second_enst, second_tpm = sorted_transcripts[1]

                # Minimum ifade kontrolü
                if dominant_tpm >= min_expression:
                    enrichment = dominant_tpm / second_tpm if second_tpm > 0 else dominant_tpm

                    # Çıktı dosyasına yaz
                    if output_handle:
                        output_handle.write(
                            f"{sample}\t{ensg}\t{dominant_enst}\t{second_enst}\t"
                            f"{dominant_tpm:.3f}\t{second_tpm:.3f}\t{enrichment:.3f}\n"
                        )

                    # Eşik kontrolü
                    if enrichment >= min_enrichment:
                        enrichment_results[ensg][dominant_enst][sample] = enrichment
                        all_enrichment_values.append(enrichment)

            processed_genes += 1
            if verbose and processed_genes % 1000 == 0:
                logger.info(f"{processed_genes} gen işlendi...")

    finally:
        if output_handle:
            output_handle.close()

    if verbose:
        logger.info(f"Zenginleşme analizi tamamlandı. Toplam {len(all_enrichment_values)} zenginleşme değeri bulundu.")

    return enrichment_results, all_enrichment_values


def read_canonical_proteins(filepath, verbose=False):
    """Canonical protein bilgilerini oku."""
    canonical_map = {}

    if verbose:
        logger.info(f"Canonical protein dosyası okunuyor: {filepath}")

    try:
        with gzip.open(filepath, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith("ENS"):
                    continue

                parts = line.split("\t")
                if len(parts) < 6:
                    continue

                canon_ensp, ensg, ensp, enst, gene_name, source = parts[:6]

                canonical_map[enst] = canon_ensp
                canonical_map[ensg] = canon_ensp
                canonical_map[ensp] = ensg

        if verbose:
            logger.info(f"Toplam {len(canonical_map)} canonical eşleşme okundu")

        return canonical_map

    except Exception as e:
        logger.error(f"Canonical dosyası okuma hatası: {e}")
        raise


def read_isoform_interactions(filepath, min_string_score, verbose=False):
    """İzoform etkileşim verilerini oku."""
    missed_interactions = defaultdict(dict)
    interaction_counts = defaultdict(int)

    if verbose:
        logger.info(f"İzoform etkileşim dosyası okunuyor: {filepath}")

    try:
        with gzip.open(filepath, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith("EN"):
                    continue

                parts = line.split("\t")
                if len(parts) < 10:
                    continue

                enst1 = parts[2]
                missed_field = parts[9]
                existing_field = parts[8]

                interaction_counts[enst1] = 0

                # Kayıp etkileşimleri işle
                for interaction_str in missed_field.split(","):
                    if ":" in interaction_str:
                        fields = interaction_str.split(":")
                        if len(fields) >= 7:
                            try:
                                string_score = float(fields[6])
                                if string_score >= min_string_score:
                                    missed_interactions[enst1][interaction_str] = 1
                                    interaction_counts[enst1] += 1
                            except (ValueError, IndexError):
                                continue

                # Mevcut etkileşimleri sayma
                for interaction_str in existing_field.split(","):
                    if ":" in interaction_str:
                        fields = interaction_str.split(":")
                        if len(fields) >= 7:
                            try:
                                string_score = float(fields[6])
                                if string_score >= min_string_score:
                                    interaction_counts[enst1] += 1
                            except (ValueError, IndexError):
                                continue

        if verbose:
            logger.info(f"Toplam {len(missed_interactions)} transkript için etkileşim verisi okundu")

        return missed_interactions, interaction_counts

    except Exception as e:
        logger.error(f"İzoform etkileşim dosyası okuma hatası: {e}")
        raise


def read_protein_sequences(filepath, verbose=False):
    """Protein sekans bilgilerini oku."""
    sequences = {}

    if verbose:
        logger.info(f"Protein sekans dosyası okunuyor: {filepath}")

    try:
        with gzip.open(filepath, "rt") as f:
            for line in f:
                line = line.strip()
                if not line or not line.startswith("EN"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 4:
                    ensp, ensg, enst, sequence = parts[:4]
                    sequences[enst] = sequence

        if verbose:
            logger.info(f"Toplam {len(sequences)} protein sekansi okundu")

        return sequences

    except Exception as e:
        logger.error(f"Protein sekans dosyası okuma hatası: {e}")
        raise


def perform_statistical_analysis(enrichment_pcawg, rel_pcawg, rel_gtex_strings,
                                 enrichment_gtex, sequences, missed_interactions,
                                 interaction_counts, canonical_map, ensg_pcawg,
                                 all_enrichment_pcawg, args):
    """Ana istatistiksel analizi gerçekleştir."""

    results = []
    pvalues = []

    # Çıktı dosyalarını aç
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

    try:
        processed_comparisons = 0

        for ensg in sorted(enrichment_pcawg):
            if ensg not in rel_gtex_strings:
                logger.debug(f"Gen GTEx'te bulunamadı: {ensg}")
                continue

            for enst in sorted(enrichment_pcawg[ensg]):
                if enst not in rel_gtex_strings[ensg]:
                    logger.debug(f"Transkript GTEx'te bulunamadı: {ensg}-{enst}")
                    continue

                for sample in sorted(enrichment_pcawg[ensg][enst]):
                    enrichment_cancer = enrichment_pcawg[ensg][enst][sample]
                    rel_expr_cancer = rel_pcawg[ensg][enst][sample]

                    # GTEx'te MDT olarak geçiyorsa atla
                    if enst in enrichment_gtex.get(ensg, {}):
                        continue

                    # GTEx'te hiç MDT yoksa atla
                    if not enrichment_gtex.get(ensg):
                        continue

                    # GTEx verilerini işle
                    gtex_expr_str = rel_gtex_strings[ensg][enst]
                    try:
                        gtex_values = [float(x) for x in gtex_expr_str.strip().split() if x.strip()]
                    except (ValueError, AttributeError):
                        logger.debug(f"GTEx değerleri okunamadı: {ensg}-{enst}")
                        continue

                    if not gtex_values:
                        continue

                    median_gtex = median(gtex_values)

                    # Binom testi için sayaçlar
                    higher = sum(1 for val in gtex_values if rel_expr_cancer > val)
                    lower = sum(1 for val in gtex_values if rel_expr_cancer < val)

                    # P-değeri hesapla
                    total_comparisons = higher + lower
                    if total_comparisons > 0:
                        pval = binomtest(higher, n=total_comparisons, p=0.5, alternative="two-sided").pvalue
                    else:
                        continue

                    # P-değeri dosyasına yaz
                    if pvalue_file:
                        pvalue_file.write(
                            f"{ensg}\t{enst}\t{sample}\t{median_gtex:.3f}\t{rel_expr_cancer:.3f}\t"
                            f"{higher}\t{lower}\t{pval:.3e}\t{gtex_expr_str}\n"
                        )

                    # Rastgele karşılaştırmalar
                    if random_pvalue_file and all_enrichment_pcawg:
                        for _ in range(100):
                            try:
                                rand_enrichment = random.choice(all_enrichment_pcawg)
                                rand_higher = sum(1 for val in gtex_values if rand_enrichment > val)
                                rand_lower = sum(1 for val in gtex_values if rand_enrichment < val)
                                rand_total = rand_higher + rand_lower
                                if rand_total > 0:
                                    rand_pval = binomtest(rand_higher, n=rand_total, p=0.5).pvalue
                                    random_pvalue_file.write(
                                        f"{ensg}\t{enst}\t{sample}\t{median_gtex:.3f}\t{rand_enrichment:.3f}\t"
                                        f"{rand_higher}\t{rand_lower}\t{rand_pval:.3e}\n"
                                    )
                            except (IndexError, ValueError):
                                continue

                    # İlk filtreleme
                    if pval >= args.maxQvalue:
                        continue

                    # GTEx MDT'lerini analiz et
                    gtex_mdt_info = analyze_gtex_mdts(ensg, enst, enrichment_gtex, sequences, args)

                    if not gtex_mdt_info['should_include']:
                        continue

                    # Etkileşim bozulmalarını analiz et
                    interaction_info = analyze_interaction_disruptions(
                        enst, ensg, sample, missed_interactions, interaction_counts,
                        canonical_map, ensg_pcawg, enrichment_gtex, args
                    )

                    # Sonuç satırını oluştur
                    result_data = {
                        'ensg': ensg,
                        'canonical_ensp': canonical_map.get(enst, '-'),
                        'enst': enst,
                        'sample': sample,
                        'gtex_mdts': gtex_mdt_info['mdt_string'],
                        'total_gtex_samples': len(gtex_values),
                        'median_gtex_enrichment': gtex_mdt_info['median_enrichment'],
                        'num_gtex_mdts': gtex_mdt_info['num_mdts'],
                        'cancer_enrichment': enrichment_cancer,
                        'median_gtex_rel_expr': median_gtex,
                        'cancer_rel_expr': rel_expr_cancer,
                        'missed_interactions_str': interaction_info['missed_str'],
                        'total_interactions': interaction_info['total_interactions'],
                        'cancer_specific_missed': interaction_info['cancer_specific'],
                        'common_missed': interaction_info['common_missed']
                    }

                    results.append(result_data)
                    pvalues.append(pval)

                    processed_comparisons += 1
                    if args.verbose and processed_comparisons % 1000 == 0:
                        logger.info(f"{processed_comparisons} karşılaştırma işlendi...")

    finally:
        if pvalue_file:
            pvalue_file.close()
        if random_pvalue_file:
            random_pvalue_file.close()

    return results, pvalues


def analyze_gtex_mdts(ensg, cancer_enst, enrichment_gtex, sequences, args):
    """GTEx MDT'lerini analiz et."""
    result = {
        'should_include': False,
        'mdt_string': '-',
        'median_enrichment': '-',
        'num_mdts': 0
    }

    if ensg not in enrichment_gtex:
        return result

    gtex_mdts = {}
    all_enrichments = []

    for gtex_enst in enrichment_gtex[ensg]:
        for sample, enrichment_val in enrichment_gtex[ensg][gtex_enst].items():
            if enrichment_val >= args.minEnrichment:
                # Farklı protein sekansı kontrolü
                if sequences.get(gtex_enst) != sequences.get(cancer_enst):
                    gtex_mdts[gtex_enst] = gtex_mdts.get(gtex_enst, 0) + 1
                    all_enrichments.append(enrichment_val)

    if not gtex_mdts:
        return result

    # %50 kuralı - en az bir MDT'nin yarıdan fazla örnekte dominant olması gerek
    total_samples = len(next(iter(enrichment_gtex[ensg].values())))

    mdt_strings = []
    for mdt_enst in sorted(gtex_mdts, key=gtex_mdts.get, reverse=True):
        count = gtex_mdts[mdt_enst]
        mdt_strings.append(f"{mdt_enst}:{count}")
        if count >= 0.5 * total_samples:
            result['should_include'] = True

    if result['should_include']:
        result['mdt_string'] = ",".join(mdt_strings)
        result['median_enrichment'] = f"{median(all_enrichments):.3f}" if all_enrichments else "-"
        result['num_mdts'] = len(all_enrichments)

    return result


def analyze_interaction_disruptions(enst, ensg, sample, missed_interactions,
                                    interaction_counts, canonical_map,
                                    ensg_pcawg, enrichment_gtex, args):
    """Protein etkileşim bozulmalarını analiz et."""
    result = {
        'missed_str': '-',
        'total_interactions': interaction_counts.get(enst, -1),
        'cancer_specific': 0,
        'common_missed': 0
    }

    if enst not in missed_interactions:
        return result

    cancer_specific_missed = []
    common_missed_count = 0

    for interaction_str in missed_interactions[enst]:
        fields = interaction_str.split(":")
        if len(fields) < 8:
            continue

        partner_ensp = fields[5]  # Etkileşimde olan diğer protein
        partner_ensg = canonical_map.get(partner_ensp)

        # Partner genin o örnekte ifade edilip edilmediğini kontrol et
        if not (partner_ensg and sample in ensg_pcawg.get(partner_ensg, {})):
            continue

        # GTEx MDT'lerinde de bu etkileşim kayboluyor mu?
        found_in_gtex = False
        if ensg in enrichment_gtex:
            for alt_enst in enrichment_gtex[ensg]:
                for alt_sample in enrichment_gtex[ensg][alt_enst]:
                    if (enrichment_gtex[ensg][alt_enst][alt_sample] >= args.minEnrichment and
                            interaction_str in missed_interactions.get(alt_enst, {})):
                        found_in_gtex = True
                        break
                if found_in_gtex:
                    break

        if found_in_gtex:
            common_missed_count += 1
        else:
            cancer_specific_missed.append(interaction_str)

    if cancer_specific_missed:
        result['missed_str'] = ",".join(cancer_specific_missed)

    result['cancer_specific'] = len(cancer_specific_missed)
    result['common_missed'] = common_missed_count

    return result


def main():
    """Ana fonksiyon."""
    args = parse_arguments()

    # Verbose modu ayarla
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("PCAWG-GTEx dominant transkript analizi başlatılıyor...")

    try:
        # Dosya kontrolü
        validate_files(args)

        # Eşleşme dosyalarını oku
        logger.info("Eşleşme dosyaları okunuyor...")
        enst2ensg = read_enst2ensg_mapping(args.ensgFile, verbose=args.verbose)
        redundant_map = read_redundant_mapping(args.redundantFile, verbose=args.verbose)
        canonical_map = read_canonical_proteins(args.canonFile, verbose=args.verbose)
        sequences = read_protein_sequences(args.sequenceFile, verbose=args.verbose)

        # Etkileşim verilerini oku
        logger.info("Protein etkileşim verileri okunuyor...")
        missed_interactions, interaction_counts = read_isoform_interactions(
            args.isoformIntFile, args.minStringScore, verbose=args.verbose
        )

        # İfade verilerini oku
        logger.info("PCAWG ifade verileri okunuyor...")
        ensg_pcawg, enst_pcawg, rel_pcawg, rel_pcawg_str, all_rel_pcawg = read_expression_data(
            args.pcawg, enst2ensg, redundant_map, args.enstColumnNuo, verbose=args.verbose
        )

        logger.info("GTEx ifade verileri okunuyor...")
        ensg_gtex, enst_gtex, rel_gtex, rel_gtex_str, all_rel_gtex = read_expression_data(
            args.gtex, enst2ensg, redundant_map, args.enstColumnNuo, verbose=args.verbose
        )

        # Zenginleşme analizleri
        logger.info("PCAWG zenginleşme analizi yapılıyor...")
        enrichment_pcawg, all_enrichment_pcawg = compute_enrichment_analysis(
            enst_pcawg, args.minExpress, args.minEnrichment,
            output_file=args.pcawgEnrichmentFile, verbose=args.verbose
        )

        logger.info("GTEx zenginleşme analizi yapılıyor...")
        enrichment_gtex, all_enrichment_gtex = compute_enrichment_analysis(
            enst_gtex, args.minExpress * 0.1, args.minEnrichment,
            output_file=args.gtexEnrichmentFile, verbose=args.verbose
        )

        # Ana istatistiksel analiz
        logger.info("İstatistiksel karşılaştırmalar yapılıyor...")
        results, pvalues = perform_statistical_analysis(
            enrichment_pcawg, rel_pcawg, rel_gtex_str, enrichment_gtex,
            sequences, missed_interactions, interaction_counts, canonical_map,
            ensg_pcawg, all_enrichment_pcawg, args
        )

        # Çoklu test düzeltmesi ve çıktı
        logger.info("Çoklu test düzeltmesi uygulanıyor...")
        if pvalues:
            _, qvalues, _, _ = multipletests(pvalues, alpha=args.maxQvalue, method='fdr_bh')

            # Header yazdır
            header = (
                "#ENSG\tENSPcanon\tcMDT\tRNAseqAliquotID\t"
                "GTExMDT\tTotalRelevantGTExSamples\tGTExMDTmedianEnrichment\tPvalue\tQvalue\t"
                "NumberOfGTExMDT\tcMDTenrichment\tGTExMDTmedianRelExpression\tcMDTrelExpression\t"
                "cMDTuniqMissedInt\tTotalNumberOfSTRINGint\t"
                "cMDTnumberOfUniqMissedInt\tNumberOfCommonMissedInt"
            )
            print(header)

            # Anlamlı sonuçları yazdır
            significant_count = 0
            for i, (result, pval, qval) in enumerate(zip(results, pvalues, qvalues)):
                if qval < args.maxQvalue:
                    output_line = (
                        f"{result['ensg']}\t{result['canonical_ensp']}\t{result['enst']}\t{result['sample']}\t"
                        f"{result['gtex_mdts']}\t{result['total_gtex_samples']}\t{result['median_gtex_enrichment']}\t"
                        f"{pval:.3e}\t{qval:.3e}\t{result['num_gtex_mdts']}\t{result['cancer_enrichment']:.3f}\t"
                        f"{result['median_gtex_rel_expr']:.3f}\t{result['cancer_rel_expr']:.3f}\t"
                        f"{result['missed_interactions_str']}\t{result['total_interactions']}\t"
                        f"{result['cancer_specific_missed']}\t{result['common_missed']}"
                    )
                    print(output_line)
                    significant_count += 1

            logger.info(
                f"Analiz tamamlandı. Toplam {len(results)} karşılaştırma, {significant_count} anlamlı sonuç bulundu.")

        else:
            logger.warning("Hiç anlamlı sonuç bulunamadı.")

    except Exception as e:
        logger.error(f"Analiz sırasında hata oluştu: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()