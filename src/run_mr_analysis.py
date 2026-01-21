import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse
import sys
import os


def run_mr_analysis(sd_asm_file, gtex_file, output_file, n_sd_default=111, n_gtex_default=500,
                    instrument_pval_threshold=0.1):
    print("----------------------------------------------------------")
    print("Starting Mendelian Randomization (Wald Ratio + Steiger Filtering)")
    print(f"  SD-ASM File: {sd_asm_file}")
    print(f"  GTEx File:   {gtex_file}")
    print(f"  Output File: {output_file}")
    print("----------------------------------------------------------")

    # 1. Check Files
    if not os.path.exists(sd_asm_file):
        print(f"Error: File not found: {sd_asm_file}")
        sys.exit(1)
    if not os.path.exists(gtex_file):
        print(f"Error: File not found: {gtex_file}")
        sys.exit(1)

    # 2. Load Data
    print("Loading SD-ASM data (Exposure)...")
    try:
        sd_asm = pd.read_csv(sd_asm_file)
        sd_asm['snp_id'] = sd_asm['snp_id'].astype(str)
        # Handle N
        if 'N' not in sd_asm.columns:
            sd_asm['N'] = n_sd_default
    except Exception as e:
        print(f"Error loading SD-ASM file: {e}")
        sys.exit(1)

    required_sd_cols = ['beta', 'varbeta', 'snp_id']
    if not all(col in sd_asm.columns for col in required_sd_cols):
        print(f"Error: SD-ASM file missing required columns: {required_sd_cols}")
        sys.exit(1)

    print("Loading GTEx data (Outcome)...")
    try:
        gtex = pd.read_csv(gtex_file)
        if 'gtex_snp_id' in gtex.columns:
            gtex['gtex_snp_id'] = gtex['gtex_snp_id'].astype(str)
        elif 'snp_id' in gtex.columns:
            gtex.rename(columns={'snp_id': 'gtex_snp_id'}, inplace=True)
            gtex['gtex_snp_id'] = gtex['gtex_snp_id'].astype(str)
        else:
            print("Error: GTEx file must contain 'gtex_snp_id' or 'snp_id' column.")
            sys.exit(1)

        # Handle N for GTEx
        if 'N' not in gtex.columns:
            gtex['N'] = n_gtex_default
        else:
            gtex['N'] = gtex['N'].fillna(n_gtex_default)

    except Exception as e:
        print(f"Error loading GTEx file: {e}")
        sys.exit(1)

    # 3. Merge Datasets
    print("Merging datasets...")
    merged_data = pd.merge(
        gtex,
        sd_asm,
        left_on='gtex_snp_id',
        right_on='snp_id',
        suffixes=('.gtex', '.sd')
    )

    print(f"Found {len(merged_data)} overlapping SNP-Gene pairs.")

    if len(merged_data) == 0:
        print("No overlapping SNPs found. Check SNP ID formats.")
        sys.exit(0)

    # 4. Filter for Valid Instruments
    merged_data['se.sd'] = np.sqrt(merged_data['varbeta.sd'])
    merged_data['z.sd'] = merged_data['beta.sd'] / merged_data['se.sd']
    merged_data['pval.sd'] = 2 * stats.norm.sf(np.abs(merged_data['z.sd']))

    valid_instruments = merged_data[merged_data['pval.sd'] < instrument_pval_threshold].copy()

    print(f"Retained {len(valid_instruments)} pairs with significant SD-ASM signal (P < {instrument_pval_threshold}).")

    if len(valid_instruments) == 0:
        print("No valid instruments found. Cannot perform MR.")
        sys.exit(0)

    # 5. Calculate Wald Ratio
    b_exp = valid_instruments['beta.sd']
    se_exp = valid_instruments['se.sd']
    b_out = valid_instruments['beta.gtex']
    se_out = np.sqrt(valid_instruments['varbeta.gtex'])

    valid_instruments['beta_causal'] = b_out / b_exp

    term1 = (se_out / b_exp) ** 2
    term2 = ((b_out * se_exp) / (b_exp ** 2)) ** 2
    valid_instruments['se_causal'] = np.sqrt(term1 + term2)

    z_mr = valid_instruments['beta_causal'] / valid_instruments['se_causal']
    valid_instruments['pval_causal'] = 2 * stats.norm.sf(np.abs(z_mr))

    # --- Steiger Filtering (Directionality Test) ---
    # Calculate Variance Explained (R2) for Exposure and Outcome
    # R2 approx = Z^2 / (N + Z^2)

    # Exposure (Methylation)
    z_exp = valid_instruments['z.sd']
    n_exp = valid_instruments['N.sd']
    valid_instruments['r2_exposure'] = (z_exp ** 2) / (n_exp + z_exp ** 2)

    # Outcome (Expression)
    z_out = b_out / se_out
    n_out = valid_instruments['N.gtex']
    valid_instruments['r2_outcome'] = (z_out ** 2) / (n_out + z_out ** 2)

    # Steiger Condition: R2_exposure > R2_outcome
    # Implies SNP -> Methylation -> Expression is more likely than SNP -> Expression -> Methylation
    valid_instruments['steiger_dir'] = valid_instruments['r2_exposure'] > valid_instruments['r2_outcome']

    # --- FDR Correction ---
    pvals = valid_instruments['pval_causal'].values
    n_tests = len(pvals)
    if n_tests > 0:
        sorted_indices = np.argsort(pvals)
        sorted_p = pvals[sorted_indices]
        ranks = np.arange(1, n_tests + 1)
        adjusted_p = sorted_p * n_tests / ranks
        adjusted_p = np.minimum.accumulate(adjusted_p[::-1])[::-1]
        adjusted_p[adjusted_p > 1] = 1.0
        fdr_values = np.zeros(n_tests)
        fdr_values[sorted_indices] = adjusted_p
        valid_instruments['fdr_causal'] = fdr_values
    else:
        valid_instruments['fdr_causal'] = []

    # 6. Organize Output
    output_cols = [
        'input_region', 'target_gene_id', 'gtex_snp_id',
        'beta.sd', 'pval.sd', 'beta.gtex',
        'beta_causal', 'se_causal', 'pval_causal', 'fdr_causal',
        'r2_exposure', 'r2_outcome', 'steiger_dir'
    ]

    rename_map = {
        'gtex_snp_id': 'snp_id',
        'beta.sd': 'beta_methylation',
        'pval.sd': 'pval_methylation',
        'beta.gtex': 'beta_expression'
    }

    final_cols = [c for c in output_cols if c in valid_instruments.columns]
    final_results = valid_instruments[final_cols].rename(columns=rename_map)
    final_results = final_results.sort_values(by='pval_causal')

    # 7. Save
    final_results.to_csv(output_file, index=False)
    print(f"\nSuccess! MR Results saved to {output_file}")

    # 8. Summary
    total_tests = len(final_results)

    # Independent Counts
    sig_fdr = len(final_results[final_results['fdr_causal'] < 0.05])
    valid_direction = len(final_results[final_results['steiger_dir'] == True])

    # Combined Count
    sig_and_valid = len(final_results[
                            (final_results['fdr_causal'] < 0.05) &
                            (final_results['steiger_dir'] == True)
                            ])

    print("\n----------------------------------------------------------")
    print("SUMMARY")
    print("----------------------------------------------------------")
    print(f"Total Associations Tested:        {total_tests}")
    print(f"Significant FDR (<0.05):          {sig_fdr}")
    print(f"Valid Causal Direction (Steiger): {valid_direction}")
    print(f"Significant + Valid Direction:    {sig_and_valid}")
    print("----------------------------------------------------------")

    print("\nTop Causal Findings (Sorted by P-value):")
    preview_cols = ['target_gene_id', 'snp_id', 'beta_causal', 'fdr_causal', 'steiger_dir']
    print(final_results[[c for c in preview_cols if c in final_results.columns]].head().to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform Mendelian Randomization (Wald Ratio + Steiger)')
    parser.add_argument('sd_asm_file', nargs='?', default='Adipocytes.sd_asm_beta.csv', help='SD-ASM CSV')
    parser.add_argument('gtex_file', nargs='?', default='Adipocytes.GTEx.summary_stats_near_sd.csv', help='GTEx CSV')
    parser.add_argument('output_file', nargs='?', default='Adipocytes_MR_Results.csv', help='Output CSV')
    parser.add_argument('--n_sd', type=int, default=111, help='Sample size for SD-ASM (Exposure)')
    parser.add_argument('--n_gtex', type=int, default=500, help='Sample size for GTEx (Outcome)')

    args = parser.parse_args()
    run_mr_analysis(args.sd_asm_file, args.gtex_file, args.output_file, args.n_sd, args.n_gtex)