import pandas as pd
import numpy as np
import statsmodels.api as sm
import scipy.stats as stats
import argparse
import gzip
import os
import sys
import csv
import heapq
import itertools

# Try to import pysam for Tabix handling
try:
    import pysam
except ImportError:
    sys.exit("CRITICAL ERROR: This script requires 'pysam' to read Tabix-indexed VCFs.\nPlease run: pip install pysam")


# ==========================================
# 1. SETUP & ARGUMENT PARSING
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(description="Check ASM vs Ancestry genome-wide using Allele-Level Regression.")
    parser.add_argument("--tissue", type=str, required=True, help="Tissue group name (e.g., 'Adipocytes', 'Endothel')")
    parser.add_argument("--output", type=str, required=True, help="Output CSV file path")

    # Paths
    parser.add_argument("--metadata", type=str, default="/cs/cbio/netanel/atlas/groups207.csv",
                        help="Path to groups metadata csv")
    parser.add_argument("--vcf_dir", type=str, default="/cs/zbio/jrosensk/atlas_genotype/vcfs",
                        help="Directory containing individual .vcf.gz files")
    parser.add_argument("--pca", type=str, default="/cs/zbio/jrosensk/atlas_genotype/vcfs/atlas_ancestry_pca.eigenvec", help="Path to PLINK .eigenvec file")
    parser.add_argument("--asm_dir", type=str, default="/cs/zbio/jrosensk/atlas_blocks_hg38_2/homog/homog_aligned/all_snps/",
                        help="Directory containing .homog.pval.gz files")
    parser.add_argument("--chart_dir", type=str,
                        default="/cs/zbio/jrosensk/atlas_blocks_hg38_2/homog/homog_aligned/all_snps/meqtl/simple_sd_asm",
                        help="Directory containing {tissue}.sd_asm_chart.1_b_bl_removed.tsv.gz files")

    return parser.parse_args()


# ==========================================
# 2. HELPER CLASSES & FUNCTIONS
# ==========================================

def asm_file_generator(filepath, sample_name):
    """
    Yields parsed lines from an ASM file.
    Yields: (chr_str, pos_int, sample_name, meth1, unmeth1, meth2, unmeth2)
    """
    if not os.path.exists(filepath):
        return

    try:
        with gzip.open(filepath, 'rt') as f:
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) < 10: continue

                # Normalize Chromosome (remove 'chr')
                c = cols[0].replace('chr', '')
                try:
                    p = int(cols[1])
                    # Parse counts: u1(4), x1(5), m1(6), u2(7), x2(8), m2(9)
                    # We need separated counts for Allele 1 and Allele 2 now
                    m1 = int(cols[6])
                    u1 = int(cols[4])

                    m2 = int(cols[9])
                    u2 = int(cols[7])

                    # Only yield if at least one allele has coverage
                    if (m1 + u1 + m2 + u2) > 0:
                        yield (c, p, sample_name, m1, u1, m2, u2)
                except ValueError:
                    continue
    except Exception as e:
        print(f"Warning: Error reading ASM file for {sample_name}: {e}")
        return


def get_genotype_from_vcf(vcf_obj, chrom, pos):
    """
    Fetch genotype from a specific individual's VCF using pysam.
    Returns: 0 (Ref/Ref), 1 (Het), 2 (Alt/Alt), or None (Missing/Error)
    """
    search_chroms = [chrom, f"chr{chrom}"] if "chr" not in chrom else [chrom, chrom.replace("chr", "")]

    rec = None
    for c in search_chroms:
        try:
            for record in vcf_obj.fetch(c, pos - 1, pos):
                if record.pos == pos:
                    rec = record
                    break
            if rec: break
        except ValueError:
            continue

    if not rec or len(rec.samples) == 0:
        return None

    gt = rec.samples[0]['GT']
    if gt is None or None in gt:
        return None

    return sum(gt)


def load_significant_loci(tissue, chart_dir):
    filename = f"{tissue}.sd_asm_chart.1_b_bl_removed.tsv.gz"
    filepath = os.path.join(chart_dir, filename)

    print(f"\n[0/4] Loading significant sites from chart: {filepath}")

    if not os.path.exists(filepath):
        sys.exit(f"ERROR: Chart file not found: {filepath}")

    try:
        df = pd.read_csv(filepath, sep='\t', header=None, compression='gzip')
        df_sig = df[df[10] < 0.05].copy()
        df_sig[0] = df_sig[0].astype(str).str.replace('chr', '', regex=False)
        sig_sites = set(zip(df_sig[0], df_sig[1]))

    except Exception as e:
        sys.exit(f"Error reading chart file with pandas: {e}")

    print(f"      Found {len(sig_sites)} significant sites (FDR < 0.05) to analyze.")
    return sig_sites


# ==========================================
# 3. MAIN LOGIC
# ==========================================
def main():
    args = parse_args()
    TISSUE = args.tissue

    print(f"--- STARTING GENOME-WIDE ANALYSIS (Allele-Level Regression) ---")
    print(f"Target Tissue: {TISSUE}")

    # --- 0. Load Target Loci ---
    target_loci = load_significant_loci(TISSUE, args.chart_dir)
    if len(target_loci) == 0:
        sys.exit("No significant sites found in chart file. Exiting.")

    # --- 1. Load Metadata ---
    print(f"\n[1/4] Loading Metadata...")
    try:
        meta_all = pd.read_csv(args.metadata)
        meta_tissue = meta_all[meta_all['group'] == TISSUE].copy()

        if len(meta_tissue) == 0:
            sys.exit(f"ERROR: No samples found for tissue group '{TISSUE}'.")

        print(f"      Found {len(meta_tissue)} samples for group '{TISSUE}'.")
        meta_tissue = meta_tissue.dropna(subset=['patientID'])
        valid_samples = meta_tissue[['name', 'patientID']].to_dict('records')

    except Exception as e:
        sys.exit(f"Error reading metadata: {e}")

    # --- 2. Load PCA ---
    print(f"\n[2/4] Loading PCA...")
    try:
        pca_df = pd.read_csv(args.pca, sep='\s+', comment='#', header=None)
        pca_map = {}
        for _, row in pca_df.iterrows():
            pid = str(row[1])
            pca_map[pid] = {'PC1': row[2], 'PC2': row[3]}
            pca_map[pid.replace('/', '.')] = {'PC1': row[2], 'PC2': row[3]}
    except Exception as e:
        sys.exit(f"Error reading PCA file: {e}")

    # --- 3. Open VCF Handles ---
    print(f"\n[3/4] Opening VCF handles...")
    vcf_handles = {}
    for sample in valid_samples:
        stub = sample['name']
        pid = str(sample['patientID'])
        vcf_path = os.path.join(args.vcf_dir, f"{pid}.vcf.gz")
        if not os.path.exists(vcf_path):
            vcf_path = os.path.join(args.vcf_dir, f"{pid.replace('/', '.')}.vcf.gz")
        if os.path.exists(vcf_path):
            try:
                vcf_handles[stub] = pysam.VariantFile(vcf_path)
            except:
                pass
    print(f"      Successfully indexed {len(vcf_handles)} VCF files.")

    # --- 4. Iterate ASM Files ---
    print(f"\n[4/4] Scanning Allele-Specific Data...")

    iterators = []
    for sample in valid_samples:
        stub = sample['name']
        if stub in vcf_handles:
            asm_path = os.path.join(args.asm_dir, f"{stub}.allele.homog.pval.gz")
            iterators.append(asm_file_generator(asm_path, stub))

    merged_stream = heapq.merge(*iterators, key=lambda x: (x[0], x[1]))
    grouped_stream = itertools.groupby(merged_stream, key=lambda x: (x[0], x[1]))

    f_out = open(args.output, 'w', newline='')
    writer = csv.writer(f_out, delimiter='\t')
    writer.writerow(
        ['Chr', 'Pos', 'N_Alleles', 'P_OLS_Ancestry_PC1', 'P_OLS_Ancestry_PC2', 'Slope_PC1', 'Slope_PC2',
         'Slope_Allele'])

    count_processed = 0
    count_hits = 0

    for (chrom, pos), group_iter in grouped_stream:
        if (chrom, pos) not in target_loci:
            for _ in group_iter: pass
            continue

        site_samples = list(group_iter)
        analysis_data = []

        for record in site_samples:
            # record: (chr, pos, stub, m1, u1, m2, u2)
            stub = record[2]
            m1, u1 = record[3], record[4]
            m2, u2 = record[5], record[6]

            # Lookup Ancestry
            patient_id = str([s['patientID'] for s in valid_samples if s['name'] == stub][0])
            if patient_id not in pca_map:
                patient_id = patient_id.replace('/', '.')
                if patient_id not in pca_map: continue
            pc1 = pca_map[patient_id]['PC1']
            pc2 = pca_map[patient_id]['PC2']

            # Lookup Genotype (to confirm phasing expectation)
            # We assume ASM file Allele 1 is Ref and Allele 2 is Alt (standard pipeline)
            # If genotype is 0/0, both are Ref. If 0/1, 1 is Ref, 2 is Alt.
            gt = get_genotype_from_vcf(vcf_handles[stub], chrom, pos)

            if gt is not None:
                # UNFOLD THE DATA: 1 Sample -> 2 Allele Observations

                # Observation 1: Allele 1 (Standard assumption: Ref-aligned)
                if (m1 + u1) > 4:
                    # Allele 1 = REF sequence. Allele 2 = ALT sequence.
                    # We code Ref as 0, Alt as 1.
                    analysis_data.append({
                        'Chr': chrom,
                        'Pos': pos,
                        'Is_Alt_Allele': 0,  # Allele 1 is Ref
                        'Meth_Ratio': m1 / (m1 + u1),
                        'PC1': pc1,
                        'PC2': pc2
                    })

                # Observation 2: Allele 2 (Standard assumption: Alt/Variant)
                if (m2 + u2) > 4:
                    analysis_data.append({
                        'Chr': chrom,
                        'Pos': pos,
                        'Is_Alt_Allele': 1,  # Allele 2 is Alt
                        'Meth_Ratio': m2 / (m2 + u2),
                        'PC1': pc1,
                        'PC2': pc2
                    })

        # --- RUN STATS (ALLELE LEVEL) ---
        # We need variance in 'Is_Alt_Allele' (meaning we need at least one Het, or one HomoRef and one HomoAlt)
        # We need variance in PCs (need more than 1 sample)

        if len(analysis_data) >= 4:  # Arbitrary minimum, 2 samples * 2 alleles
            df = pd.DataFrame(analysis_data)

            # Check for variance
            if df['Meth_Ratio'].std() > 0 and df['Is_Alt_Allele'].std() > 0:
                try:
                    # Standardize for Beta coefficients
                    df_std = df[['PC1', 'PC2', 'Is_Alt_Allele']].copy()

                    for col in df_std.columns:
                        if df_std[col].std() > 0:
                            df_std[col] = (df_std[col] - df_std[col].mean()) / df_std[col].std()
                        else:
                            df_std[col] = 0.0

                    X = df_std
                    X = sm.add_constant(X)
                    Y = df['Meth_Ratio']

                    model = sm.OLS(Y, X).fit()
                    p_PC1 = model.pvalues.get('PC1', "NA")
                    p_PC2 = model.pvalues.get('PC2', "NA")
                    slope_PC1 = model.params.get('PC1', "NA")
                    slope_PC2 = model.params.get('PC2', "NA")
                    slope_GT = model.params.get('Is_Alt_Allele', "NA")

                    writer.writerow([
                        chrom, pos, len(df), p_PC1, p_PC2, slope_PC1, slope_PC2, slope_GT
                    ])
                    count_hits += 1
                except:
                    pass
            else:
                # No variance (e.g. all alleles 0, or all Meth 100%)
                writer.writerow([chrom, pos, len(df), "NA", "NA", "NA", "NA", "NA"])

        count_processed += 1
        if count_processed % 100 == 0:
            print(f"\rProcessed {count_processed} significant sites... (Hits: {count_hits})", end="")

    print(f"\n[5/5] Done. Processed {count_processed} significant sites. Results: {args.output}")
    f_out.close()

    # --- Summary Report ---
    print(f"\n[6/5] Generating Summary Report...")
    try:
        results = pd.read_csv(args.output, sep='\t')
        results['Slope_Allele_Std'] = pd.to_numeric(results['Slope_Allele'], errors='coerce')
        valid_results = results.dropna(subset=['Slope_Allele_Std']).copy()

        # Convert columns
        cols = ['P_OLS_Ancestry_PC1', 'P_OLS_Ancestry_PC2', 'Slope_PC1', 'Slope_PC2']
        for c in cols: valid_results[c] = pd.to_numeric(valid_results[c], errors='coerce')

        N = len(valid_results)
        if N > 0:
            sig_pc1 = len(valid_results[valid_results['P_OLS_Ancestry_PC1'] < 0.05])
            sig_pc2 = len(valid_results[valid_results['P_OLS_Ancestry_PC2'] < 0.05])

            # Compare Absolute Standardized Slopes
            pc1_bigger = len(
                valid_results[valid_results['Slope_PC1'].abs() > valid_results['Slope_Allele'].abs()])
            allele_bigger_than_pc1 = len(
                valid_results[valid_results['Slope_Allele'].abs() > valid_results['Slope_PC1'].abs()])

            pc2_bigger = len(
                valid_results[valid_results['Slope_PC2'].abs() > valid_results['Slope_Allele'].abs()])
            allele_bigger_than_pc2 = len(
                valid_results[valid_results['Slope_Allele'].abs() > valid_results['Slope_PC2'].abs()])

            print("\n" + "=" * 60)
            print(f"SUMMARY STATISTICS (N={N} analyzed sites)")
            print("Technique: Allele-Level Regression (Unfolded)")
            print("Slope 'Allele' represents effect of sequence (Ref vs Alt)")
            print("=" * 60)
            print(f"Ancestry Significance (P < 0.05):")
            print(f"  PC1 Significant: {sig_pc1} ({sig_pc1 / N * 100:.1f}%)")
            print(f"  PC2 Significant: {sig_pc2} ({sig_pc2 / N * 100:.1f}%)")
            print("-" * 60)
            print(f"Effect Size Comparison (Absolute Standardized Slopes):")
            print(f"  |Slope PC1| > |Slope Allele|: {pc1_bigger} ({pc1_bigger / N * 100:.1f}%)")
            print(f"  |Slope Allele| > |Slope PC1|: {allele_bigger_than_pc1} ({allele_bigger_than_pc1 / N * 100:.1f}%)")
            print("-" * 60)
            print(f"  |Slope PC2| > |Slope Allele|: {pc2_bigger} ({pc2_bigger / N * 100:.1f}%)")
            print(f"  |Slope Allele| > |Slope PC2|: {allele_bigger_than_pc2} ({allele_bigger_than_pc2 / N * 100:.1f}%)")
            print("=" * 60)
    except Exception as e:
        print(f"Error summary: {e}")


if __name__ == "__main__":
    main()