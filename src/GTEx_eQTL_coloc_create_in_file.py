import pandas as pd
import numpy as np
import argparse
import sys
import os


def process_chunk(df):
    """
    Helper function to process a single chunk of GTEx data.
    """
    # 1. Check for required columns
    # We allow flexible column names if 'pval_nominalslope' appears due to header issues
    if 'slope' not in df.columns:
        # Fallback for weird headers or specific GTEx versions
        print("Warning: 'slope' column not found. Checking for alternative headers...")
        if 'pval_nominalslope' in df.columns:
            print(
                "Error: Your file header seems corrupted (merged columns). Please ensure it is a valid tab-separated file.")
            sys.exit(1)
        # If standard columns are missing
        print(f"Available columns: {list(df.columns)}")
        raise ValueError("Could not find 'slope' or 'slope_se' columns.")

    # 2. Extract Beta (Slope) and Standard Error
    out_df = pd.DataFrame()
    out_df['beta'] = df['slope']
    out_df['se'] = df['slope_se']
    out_df['varbeta'] = df['slope_se'] ** 2

    # 3. Format SNP ID
    # GTEx format: chr1_12345_A_G_b38
    # Target format: chr1:12345
    # We split the string by '_' and take the first two parts
    # Vectorized string split is faster
    split_ids = df['variant_id'].astype(str).str.split('_', expand=True)
    out_df['snp_id'] = split_ids[0] + ':' + split_ids[1]

    # 4. Extract other useful info
    if 'pval_nominal' in df.columns:
        out_df['pval'] = df['pval_nominal']
    elif 'pval_beta' in df.columns:
        out_df['pval'] = df['pval_beta']

    out_df['maf'] = df['maf']
    out_df['gene_id'] = df['gene_id']  # Keep gene ID to filter later in R

    # 5. Calculate Sample Size (N) per row
    # N approx = minor_allele_count / (2 * maf)
    # Handle division by zero or NA
    out_df['N'] = df['minor_allele_count'] / (2 * df['maf'])
    out_df['N'] = out_df['N'].round().fillna(0).astype(int)

    return out_df


def process_gtex_file(input_file, output_file, chunk_size=100000):
    print(f"Reading file in chunks of {chunk_size}: {input_file}")

    if os.path.exists(output_file):
        os.remove(output_file)

    chunk_iter = 0
    total_rows = 0

    try:
        # Use sep='\t' for GTEx files
        with pd.read_csv(input_file, sep='\t', chunksize=chunk_size, compression='infer') as reader:
            for chunk in reader:
                try:
                    processed_chunk = process_chunk(chunk)

                    mode = 'w' if chunk_iter == 0 else 'a'
                    header = True if chunk_iter == 0 else False

                    processed_chunk.to_csv(output_file, index=False, mode=mode, header=header)

                    total_rows += len(processed_chunk)
                    chunk_iter += 1
                    print(f"Processed chunk {chunk_iter} ({total_rows} rows total)...", end='\r')

                except ValueError as ve:
                    print(f"\nError in chunk {chunk_iter}: {ve}")
                    sys.exit(1)

        print(f"\nSuccess! Processed {total_rows} rows. Saved to {output_file}")

        # WARNING for eGenes file
        if total_rows < 50000:  # Heuristic: eGenes files are small (~20k-30k rows)
            print("\n*** WARNING ***")
            print("The output file seems small (~30k rows).")
            print("If you used the 'egenes.txt' file, you only have the TOP SNP per gene.")
            print("For Coloc, you usually need the full summary statistics ('all_pairs.txt').")
            print("***************\n")

    except Exception as e:
        print(f"\nError processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Format GTEx Summary Stats for Coloc.')
    parser.add_argument('input_file', help='Path to GTEx .txt.gz file')
    parser.add_argument('--output', '-o', default='gtex_coloc_input.csv',
                        help='Path to output CSV')

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)

    process_gtex_file(args.input_file, args.output)