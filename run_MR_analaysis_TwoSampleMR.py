import argparse
import pandas as pd
import numpy as np
import sys
import os

# RPy2 Imports
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.rinterface_lib.embedded import RRuntimeError


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run MR: Methylation (Exposure) vs GTEx Expression (Outcome)")
    parser.add_argument("--asm", required=True, help="Path to SD-ASM (Methylation) summary statistics CSV")
    parser.add_argument("--gtex", required=True, help="Path to GTEx summary statistics CSV")
    parser.add_argument("--out", required=True, help="Path for the output CSV file")
    parser.add_argument("--p_threshold", type=float, default=0.05,
                        help="P-value threshold for mQTL instrument selection (default: 0.05)")
    return parser.parse_args()


def get_position_clusters(snp_list, window_size=100):
    chr_groups = {}
    for snp in snp_list:
        try:
            parts = snp.split(':')
            pos = int(parts[-1])
            chrom = ":".join(parts[:-1])
            if chrom not in chr_groups:
                chr_groups[chrom] = []
            chr_groups[chrom].append({'id': snp, 'pos': pos})
        except ValueError:
            continue

    final_clusters = []
    for chrom, snps in chr_groups.items():
        snps.sort(key=lambda x: x['pos'])
        n, i = len(snps), 0
        while i < n:
            seed_snp = snps[i]
            current_cluster = [seed_snp['id']]
            j = i + 1
            while j < n:
                if (snps[j]['pos'] - seed_snp['pos']) <= window_size:
                    current_cluster.append(snps[j]['id'])
                    j += 1
                else:
                    break
            final_clusters.append(current_cluster)
            i = j
    return final_clusters


def run_mr_analysis(asm_file, gtex_file, output_file, p_thresh):
    print(f"Loading libraries and data...")
    try:
        tsmr = importr('TwoSampleMR')
    except Exception:
        print("Error: R package 'TwoSampleMR' not found.")
        sys.exit(1)

    try:
        asm_df = pd.read_csv(asm_file)
        gtex_df = pd.read_csv(gtex_file)
    except Exception as e:
        print(f"Error reading input: {e}");
        sys.exit(1)

    # 1. Preprocess Methylation (Now the Exposure)
    if 'varbeta' in asm_df.columns:
        asm_df['se'] = np.sqrt(asm_df['varbeta'])
    asm_df['alt'], asm_df['ref'] = 'A', 'G'

    # Standardize Methylation stats
    asm_beta_sd = asm_df['beta'].std()
    asm_df['beta_std'] = asm_df['beta'] / asm_beta_sd
    asm_df['se_std'] = asm_df['se'] / asm_beta_sd

    # 2. Preprocess GTEx (Now the Outcome)
    if 'gtex_snp_id' in gtex_df.columns:
        gtex_df.rename(columns={'gtex_snp_id': 'snp_id'}, inplace=True)
    if 'varbeta' in gtex_df.columns:
        gtex_df['se'] = np.sqrt(gtex_df['varbeta'])
    gtex_df['alt'], gtex_df['ref'] = 'A', 'G'

    # Standardize GTEx stats
    gtex_beta_sd = gtex_df['beta'].std()
    gtex_df['beta_std'] = gtex_df['beta'] / gtex_beta_sd
    gtex_df['se_std'] = gtex_df['se'] / gtex_beta_sd

    # 3. Filter for common SNPs and strong mQTLs
    common_snps = set(gtex_df['snp_id']).intersection(set(asm_df['snp_id']))
    asm_df = asm_df[asm_df['snp_id'].isin(common_snps)].copy()

    # EXPOSURE FILTERING: Keep only SNPs that are significant mQTLs
    asm_df = asm_df[asm_df['meqtl_p_val'] < p_thresh]

    if asm_df.empty:
        print(f"No SNPs passed the mQTL threshold of {p_thresh}.")
        return

    gtex_df = gtex_df[gtex_df['snp_id'].isin(asm_df['snp_id'])].copy()
    unique_genes = gtex_df['target_gene_id'].unique()
    all_results = []

    print(f"Processing {len(unique_genes)} genes with valid mQTL instruments...")

    for i, gene in enumerate(unique_genes):
        if i % 50 == 0:
            print(f"   ... gene {i + 1}/{len(unique_genes)}: {gene}")

        gene_outcome_full = gtex_df[gtex_df['target_gene_id'] == gene]
        matched_snps = gene_outcome_full['snp_id'].tolist()
        if not matched_snps: continue

        snp_clusters = get_position_clusters(matched_snps, window_size=100)

        for cluster_idx, cluster_snps in enumerate(snp_clusters):
            # Isolate data for this specific gene-cluster pair
            exposure_subset = asm_df[asm_df['snp_id'].isin(cluster_snps)].copy()
            outcome_subset = gene_outcome_full[gene_outcome_full['snp_id'].isin(cluster_snps)].copy()

            cluster_id = f"{gene}_C{cluster_idx + 1}"
            exposure_subset['pheno'], outcome_subset['pheno'] = cluster_id, 'Expression'

            try:
                with localconverter(ro.default_converter + pandas2ri.converter):
                    # Format Exposure (Methylation)
                    exp_dat = tsmr.format_data(
                        dat=exposure_subset, type='exposure', snp_col='snp_id',
                        phenotype_col='pheno', beta_col='beta_std', se_col='se_std',
                        pval_col='meqtl_p_val', effect_allele_col='alt', other_allele_col='ref'
                    )
                    # Format Outcome (Expression)
                    out_dat = tsmr.format_data(
                        dat=outcome_subset, type='outcome', snp_col='snp_id',
                        phenotype_col='pheno', beta_col='beta_std', se_col='se_std',
                        pval_col='pval', effect_allele_col='alt', other_allele_col='ref'
                    )

                    dat = tsmr.harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat, action=3)
                    res = tsmr.mr_singlesnp(dat)
                    res_py = ro.conversion.rpy2py(res)

                if not res_py.empty:
                    for _, row in res_py.iterrows():
                        if pd.isna(row.get('b')): continue

                        raw_label = str(row['SNP'])
                        if "All - Inverse variance weighted" in raw_label:
                            method, snp_display = 'IVW (Cluster)', ";".join(cluster_snps)
                        elif "All - MR Egger" in raw_label:
                            method, snp_display = 'MR Egger', ";".join(cluster_snps)
                        else:
                            method, snp_display = 'Wald Ratio', raw_label

                        all_results.append({
                            'gene_id': gene, 'cluster_id': cluster_id, 'snp_id': snp_display,
                            'mr_beta': row['b'], 'mr_se': row['se'], 'mr_pval': row['p'],
                            'method': method
                        })

            except Exception as e:
                print(f"ERROR: {e}")
                continue

    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df.to_csv(output_file, index=False)
        print(f"Done. Results written to {output_file}")
    else:
        print("No significant MR associations found.")


if __name__ == "__main__":
    args = parse_args()
    run_mr_analysis(args.asm, args.gtex, args.out, args.p_threshold)