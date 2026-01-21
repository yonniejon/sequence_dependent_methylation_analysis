import pandas as pd
import numpy as np
import argparse
import sys
import os
import scipy.stats  # Required for P-value calculation

def process_chunk(df):
    """
    Helper function to process a single chunk of data.
    Returns the processed DataFrame with beta, varbeta, calculated p-value, and original fisher p-value.
    """
    # 1. Clean header: remove '#' from '#chrom' if present
    if '#chrom' in df.columns:
        df.rename(columns={'#chrom': 'chrom'}, inplace=True)
    
    # Check for required columns in this chunk
    required_cols = ['m1', 'u1', 'm2', 'u2', 'chrom', 'start']
    # We only check columns on the first chunk usually, but good to be safe
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        # If columns are missing, we can't process. 
        # But we return None to signal failure in the main loop.
        raise ValueError(f"Missing required columns: {missing_cols}")

    # 2. Define Counts
    m1 = df['m1']
    u1 = df['u1']
    m2 = df['m2']
    u2 = df['u2']
    total = m1 + m2 + u1 + u2

    # 3. Add Pseudocounts (Haldane-Anscombe correction)
    pseudo = 0.5
    m1_adj = m1 + pseudo
    u1_adj = u1 + pseudo
    m2_adj = m2 + pseudo
    u2_adj = u2 + pseudo

    # 4. Calculate Odds Ratio (OR)
    odds_ratio = (m1_adj * u2_adj) / (u1_adj * m2_adj)

    # 5. Calculate Beta (Log Odds Ratio)
    # Create a copy to avoid SettingWithCopyWarning if df is a slice
    out_df = pd.DataFrame()
    out_df['snp_id'] = df['chrom'].astype(str) + ':' + df['start'].astype(str)
    out_df['beta'] = np.log(odds_ratio)

    # 6. Calculate Variance of Beta (Woolf's formula)
    out_df['varbeta'] = (1/m1_adj) + (1/u1_adj) + (1/m2_adj) + (1/u2_adj)

    # 7. Calculate P-value (Wald Test)
    # We calculate this from the Beta and VarBeta to ensure consistency for Coloc.
    # Z = Beta / SE
    se = np.sqrt(out_df['varbeta'])
    z_score = out_df['beta'] / se
    
    # Calculate 2-sided p-value using Survival Function (sf), which is 1 - CDF.
    # This is often more precise for very small p-values.
    out_df['meqtl_p_val'] = 2 * scipy.stats.norm.sf(np.abs(z_score))

    # 8. Include original Fisher's p-value
    # Assuming the input column is named 'meqtl_p_val'
    out_df['fisher_p-value'] = df['meqtl_p_val']
    out_df["N"] = total
        
    return out_df

def calculate_summary_stats_chunked(input_file, output_file, chunk_size=100000):
    """
    Reads SD-ASM count data in chunks to save memory.
    """
    print(f"Reading file in chunks of {chunk_size}: {input_file}")
    
    # Clear output file if it exists so we can append cleanly
    if os.path.exists(output_file):
        os.remove(output_file)
        
    chunk_iter = 0
    total_rows = 0
    
    try:
        # We use sep='\t' assuming it is a TSV. 
        # If it is variable whitespace, use delim_whitespace=True instead.
        # compression='infer' will automatically handle .gz files.
        with pd.read_csv(input_file, sep='\t', chunksize=chunk_size, compression='infer') as reader:
            for chunk in reader:
                try:
                    processed_chunk = process_chunk(chunk)
                    
                    # Write to CSV
                    # If it's the first chunk, write header. Otherwise, skip header.
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

    except Exception as e:
        print(f"\nError processing file: {e}")
        print("Tip: If the file is not tab-separated, try changing the script to use sep=r'\\s+' (note: this is slower).")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process SD-ASM files for Coloc analysis (Low Memory Mode).')
    parser.add_argument('input_file', help='Path to the input SD-ASM text file (can be .gz)')
    parser.add_argument('--output', '-o', default='sd_asm_coloc_input.csv', 
                        help='Path to the output CSV file')
    parser.add_argument('--chunksize', type=int, default=100000, 
                        help='Number of rows to process at a time (default: 100,000)')

    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
        
    calculate_summary_stats_chunked(args.input_file, args.output, args.chunksize)
