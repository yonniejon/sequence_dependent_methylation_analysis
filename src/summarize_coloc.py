import pandas as pd
import argparse
import sys
import os

def summarize_coloc_results(input_file, h4_threshold=0.5, h3_threshold=0.5):
    """
    Reads a Coloc results CSV and prints a summary of significant findings.
    """
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)

    # basic validation
    required_cols = ['input_region', 'target_gene', 'PP.H3', 'PP.H4']
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Input CSV missing required columns. Found: {list(df.columns)}")
        sys.exit(1)

    print(f"\n--- Coloc Results Summary: {input_file} ---")
    print(f"Total Regions Tested: {len(df)}")

    # 1. Strong Colocalization (H4)
    # Evidence that SD-ASM and Gene Expression share a causal variant
    coloc_hits = df[df['PP.H4'] >= h4_threshold].copy()
    coloc_hits = coloc_hits.sort_values(by='PP.H4', ascending=False)

    print(f"\n[1] Evidence for Colocalization (PP.H4 >= {h4_threshold}): {len(coloc_hits)} hits")
    if not coloc_hits.empty:
        # Create a cleaner view
        view = coloc_hits[['input_region', 'target_gene', 'n_snps', 'PP.H3', 'PP.H4']]
        print(view.to_string(index=False))
    else:
        print("  No hits found matching criteria.")

    # 2. Distinct Signals (H3)
    # Evidence that both are genetic, but driven by DIFFERENT variants
    distinct_hits = df[df['PP.H3'] >= h3_threshold].copy()
    distinct_hits = distinct_hits.sort_values(by='PP.H3', ascending=False)

    print(f"\n[2] Evidence for Distinct Signals (PP.H3 >= {h3_threshold}): {len(distinct_hits)} hits")
    if not distinct_hits.empty:
        view = distinct_hits[['input_region', 'target_gene', 'n_snps', 'PP.H3', 'PP.H4']]
        print(view.to_string(index=False))
    else:
        print("  No hits found matching criteria.")

    # 3. Top H2 (GTEx driven, SD-ASM quiet)
    # Often happens if SD-ASM lacks power or is truly not genetic
    top_h2 = df[df['PP.H2'] > 0.8]
    print(f"\n[3] Strong GTEx signal only (PP.H2 > 0.8): {len(top_h2)} regions")

    # 4. Summary Stats
    print("\n--- Statistics ---")
    print(f"Max PP.H4 observed: {df['PP.H4'].max():.4f}")
    print(f"Max PP.H3 observed: {df['PP.H3'].max():.4f}")
    print(f"Mean SNP overlap:   {df['n_snps'].mean():.1f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarize Coloc Analysis Results')
    parser.add_argument('input_file', help='Path to Coloc results CSV')
    parser.add_argument('--h4', type=float, default=0.5, help='Threshold for PP.H4 (Colocalization) [Default: 0.5]')
    parser.add_argument('--h3', type=float, default=0.5, help='Threshold for PP.H3 (Distinct Signals) [Default: 0.5]')

    args = parser.parse_args()
    summarize_coloc_results(args.input_file, args.h4, args.h3)