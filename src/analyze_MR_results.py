import pandas as pd
import argparse
from statsmodels.stats.multitest import multipletests

def analyze_mr_results(input_file, p_threshold=0.05):
    # 1. Load data
    df = pd.read_csv(input_file)
    
    # 2. Drop rows that are NA (no valid MR test performed)
    df = df.dropna(subset=['mr_beta', 'mr_pval'])
    
    # 3. Identify unique clusters and select the primary test
    # If a cluster has both Wald Ratio and IVW, we want the IVW.
    # We sort by method so 'IVW (Cluster)' comes before 'Wald Ratio'
    df['method_priority'] = df['method'].map({'IVW (Cluster)': 0, 'Wald Ratio': 1, 'MR Egger': 2})
    
    # For each unique cluster, keep only the most reliable method
    primary_tests = df.sort_values(['cluster_id', 'method_priority']).drop_duplicates('cluster_id')
    
    total_clusters = len(primary_tests)
    
    if total_clusters == 0:
        print("No valid (non-NA) clusters found in the file.")
        return

    # 4. Count nominally significant results
    nominal_sig = (primary_tests['mr_pval'] < p_threshold).sum()
    
    # 5. Apply FDR (Benjamini-Hochberg) correction
    # returns: (rejected_null, pvals_corrected, alphacSidaka, alphacBonf)
    rejected, pvals_fdr, _, _ = multipletests(primary_tests['mr_pval'], 
                                              alpha=p_threshold, 
                                              method='fdr_bh')
    
    fdr_sig = rejected.sum()
    
    # Output Results
    print("-" * 40)
    print(f"File: {input_file}")
    print("-" * 40)
    print(f"Total valid clusters tested:    {total_clusters}")
    print(f"Nominally significant (p < {p_threshold}): {nominal_sig}")
    print(f"FDR significant (q < {p_threshold}):       {fdr_sig}")
    print("-" * 40)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize MR Cluster results.")
    parser.add_argument("file", help="Path to the MR results CSV")
    parser.add_argument("--p", type=float, default=0.05, help="P-value threshold (default: 0.05)")
    
    args = parser.parse_args()
    analyze_mr_results(args.file, args.p)
