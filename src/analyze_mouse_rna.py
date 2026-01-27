import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests
from pybedtools import BedTool
import seaborn as sns
import matplotlib.pyplot as plt
import os
from sklearn.neighbors import NearestNeighbors

# --- Configuration ---
GENE_COUNTS_FILE = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/GSE106208_counts_per_gene.tpm.txt.gz'
GENES_BED_MM10 = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/genes_mm10.bed'  # The output from Part 1
DMR_FILE = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/B6_C3.DMRS.tsv'
BASE_DIR = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/'
DE_INFO_FILE = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/genes_with_DE_info.tsv'
ENHANCER_FILE = os.path.join(BASE_DIR, 'enhancer_candidates.tsv')
SILENCER_FILE = os.path.join(BASE_DIR, 'silencer_candidates.tsv')
ENHANCER_SILENCER_FILE = os.path.join(BASE_DIR, 'combined_enhancer_silencer_candidates.tsv')
OUTPUT_PLOT = os.path.join(BASE_DIR, 'B6_C3_Meth_vs_Expr_Plot.png')
WINDOW_SIZE = 150000  # 150kb

def run_matched_enrichment(de_df, de_genes_set, near_dmr_set, k=10):
    """
    Performs enrichment analysis using an expression-matched background.
    """
    print(f"\nCalculating Matched Enrichment (K={k})...")

    # 1. Prepare Expression Data
    # Calculate global average TPM for matching
    de_df['avg_tpm'] = de_df[['B6_mean', 'C3_mean']].mean(axis=1)

    # Split into DE and Non-DE pools
    de_pool = de_df[de_df['geneID'].isin(de_genes_set)].copy()
    non_de_pool = de_df[~de_df['geneID'].isin(de_genes_set)].copy()

    if de_pool.empty or non_de_pool.empty:
        print("Error: One of the gene pools is empty.")
        return

    # 2. Find Nearest Neighbors based on TPM
    # We reshape to 2D for sklearn: [[val1], [val2]...]
    non_de_tpm = non_de_pool['avg_tpm'].values.reshape(-1, 1)
    de_tpm = de_pool['avg_tpm'].values.reshape(-1, 1)

    nn = NearestNeighbors(n_neighbors=k, algorithm='ball_tree')
    nn.fit(non_de_tpm)

    # distances, indices are the positions in the non_de_pool
    distances, indices = nn.kneighbors(de_tpm)

    # Flatten indices and get unique control genes
    matched_indices = np.unique(indices.flatten())
    matched_control_genes = set(non_de_pool.iloc[matched_indices]['geneID'])

    # 3. Construct Contingency Table
    # Row 1: DE Genes
    # Row 2: Matched Non-DE Genes (Control)

    a = len(de_genes_set.intersection(near_dmr_set))
    b = len(de_genes_set - near_dmr_set)

    c = len(matched_control_genes.intersection(near_dmr_set))
    d = len(matched_control_genes - near_dmr_set)

    # 4. Fisher's Exact Test
    oddsratio, pvalue = fisher_exact([[a, b], [c, d]])

    print(f"--- Matched Enrichment Results ---")
    print(f"DE Genes (Total {len(de_genes_set)}): {a} near DMR, {b} far")
    print(f"Matched Background (Total {len(matched_control_genes)}): {c} near DMR, {d} far")
    print(f"Matched Odds Ratio: {oddsratio:.4f}")
    print(f"Fisher Exact p-value: {pvalue:.4e}")

    return oddsratio, pvalue


def enhancer_analysis(df_dmr, genes_bed, significant_genes_df):
    # Re-run window to get specific gene-DMR pairs with metadata
    # The 'genes_bed' has standard bed columns + name (geneID)
    # The 'dmr_bed' has standard bed columns + region + b6_meth + c3_meth + pval + status

    # window returns: [gene_chrom, gene_start, gene_end, gene_name, gene_score, gene_strand,
    #                  dmr_chrom, dmr_start, dmr_end, dmr_region, b6_meth, c3_meth, dmr_pval, dmr_status]
    dmr_cols = ['#chrom', 'start', 'end', 'region', 'B6 mean meth.', 'C3 mean meth.', 'p-value', 'B6 status']
    dmr_bed = BedTool.from_dataframe(df_dmr[dmr_cols])
    pairs = genes_bed.window(dmr_bed, w=WINDOW_SIZE)

    # Convert pairs to dataframe
    # Note: BedTool.window results are generic; we map indices manually based on input structure
    # Gene BED had 6 cols. DMR BED had 8 cols. Total 14.
    pairs_df = pairs.to_dataframe(names=[
        'g_chrom', 'g_start', 'g_end', 'geneID', 'g_score',
        'd_chrom', 'd_start', 'd_end', 'dmr_region', 'dmr_b6_meth', 'dmr_c3_meth', 'dmr_pval', 'dmr_status'
    ])

    # Merge with Differential Expression Results (Inner join keeps only genes that are DE)
    merged = pd.merge(pairs_df, significant_genes_df, on='geneID', how='inner')

    status_counts = merged.groupby('geneID')['dmr_status'].nunique()
    consistent_genes = status_counts[status_counts == 1].index

    merged_consistent = merged[merged['geneID'].isin(consistent_genes)].copy()

    enhancer_list = []
    silencer_list = []

    for idx, row in merged_consistent.iterrows():
        b6_med = row['B6_median']
        c3_med = row['C3_median']
        dmr_stat = row['dmr_status']  # hypermethyalted or hypomethylated

        # Classification Logic
        is_enhancer = False
        is_silencer = False

        # Enhancer: (Hyper B6 & High C3) OR (Hypo B6 & High B6) -> Negative Correlation
        if dmr_stat == 'hypermethyalted' and c3_med > b6_med:
            is_enhancer = True
        elif dmr_stat == 'hypomethylated' and b6_med > c3_med:
            is_enhancer = True

        # Silencer: (Hyper B6 & High B6) OR (Hypo B6 & High C3) -> Positive Correlation
        if dmr_stat == 'hypermethyalted' and b6_med > c3_med:
            is_silencer = True
        elif dmr_stat == 'hypomethylated' and c3_med > b6_med:
            is_silencer = True

        # Prepare output row
        out_row = {
            'geneID': row['geneID'],
            'gene_chrom': row['g_chrom'],
            'gene_start': row['g_start'],
            'gene_end': row['g_end'],
            'B6_median_expr': b6_med,
            'C3_median_expr': c3_med,
            'DE_FDR': row['fdr'],
            'DMR_region': row['dmr_region'],
            'DMR_B6_meth': row['dmr_b6_meth'],
            'DMR_C3_meth': row['dmr_c3_meth'],
            'DMR_pval': row['dmr_pval'],
            'DMR_status': dmr_stat
        }

        if is_enhancer:
            enhancer_list.append(out_row)
        if is_silencer:
            silencer_list.append(out_row)

    # Save Files
    pd.DataFrame(enhancer_list).to_csv(ENHANCER_FILE, sep='\t', index=False)
    pd.DataFrame(silencer_list).to_csv(SILENCER_FILE, sep='\t', index=False)

    df_enh = pd.DataFrame(enhancer_list)
    df_sil = pd.DataFrame(silencer_list)

    # 2. Add the labeling column
    df_enh['regulatory_class'] = 'enhancer'
    df_sil['regulatory_class'] = 'silencer'

    # 3. Combine them into one
    # ignore_index=True prevents duplicate index values from the two original DFs
    combined_df = pd.concat([df_enh, df_sil], ignore_index=True)

    # 4. Save to a single file
    combined_df.to_csv(ENHANCER_SILENCER_FILE, sep='\t', index=False)

    print(f"Enhancer candidates found: {len(enhancer_list)}")
    print(f"Silencer candidates found: {len(silencer_list)}")


def load_and_analyze():
    # ---------------------------------------------------------
    # 1. DMR Processing and Spatial Association (Same as before)
    # ---------------------------------------------------------
    print("Loading DMRs and Genes...")
    genes_bed = BedTool(GENES_BED_MM10)
    gene_coords_df = pd.read_csv(GENES_BED_MM10, sep='\t', header=None,
                                 names=['chrom', 'start', 'end', 'geneID', 'score'])
    df_dmr = pd.read_csv(DMR_FILE, sep='\t')

    hyper_df = df_dmr[df_dmr['B6 status'] == 'hypermethylated']
    hypo_df = df_dmr[df_dmr['B6 status'] == 'hypomethylated']

    hyper_bed = BedTool.from_dataframe(hyper_df[['#chrom', 'start', 'end', 'region']])
    hypo_bed = BedTool.from_dataframe(hypo_df[['#chrom', 'start', 'end', 'region']])

    # Identify genes within 150kb
    genes_near_hyper = genes_bed.window(hyper_bed, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique()
    genes_near_hypo = genes_bed.window(hypo_bed, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique()

    # Save lists
    pd.Series(genes_near_hyper).to_csv('/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/genes_near_hyper_dmrs.csv', index=False, header=['geneID'])
    pd.Series(genes_near_hypo).to_csv('/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/genes_near_hypo_dmrs.csv', index=False, header=['geneID'])

    set_hyper = set(genes_near_hyper)
    set_hypo = set(genes_near_hypo)

    common = set_hyper.intersection(set_hypo)
    unique_hyper = set_hyper - set_hypo
    unique_hypo = set_hypo - set_hyper

    total_associated = len(set_hyper.union(set_hypo))

    if total_associated > 0:
        print(f"\n--- Association Stats ---")
        print(f"Common genes: {len(common)} ({len(common) / total_associated:.2%})")
        print(f"Unique to Hyper: {len(unique_hyper)} ({len(unique_hyper) / total_associated:.2%})")
        print(f"Unique to Hypo: {len(unique_hypo)} ({len(unique_hypo) / total_associated:.2%})")
    # ---------------------------------------------------------
    # 2. Differential Expression Analysis (Updated)
    # ---------------------------------------------------------
    print("\nRunning Differential Expression Analysis...")
    if not os.path.exists(DE_INFO_FILE):
        df_counts = pd.read_csv(GENE_COUNTS_FILE, sep='\t')

        # Identify Strain Columns
        b6_cols = [c for c in df_counts.columns if 'B6' in c]
        c3_cols = [c for c in df_counts.columns if 'C3' in c]

        de_results = []

        # Iterate over genes
        for i in range(len(df_counts)):
            row = df_counts.iloc[i]
            gene_id = row['geneID']

            # Get values and ensure float type
            b6_vals = row[b6_cols].values.astype(float)
            c3_vals = row[c3_cols].values.astype(float)

            # Calculate Statistics
            b6_mean = np.mean(b6_vals)
            b6_median = np.median(b6_vals)
            c3_mean = np.mean(c3_vals)
            c3_median = np.median(c3_vals)

            # Mann-Whitney U Test
            try:
                stat, pval = mannwhitneyu(b6_vals, c3_vals, alternative='two-sided')
            except ValueError:
                pval = 1.0

            de_results.append({
                'geneID': gene_id,
                'p_value': pval,
                'B6_mean': b6_mean,
                'B6_median': b6_median,
                'C3_mean': c3_mean,
                'C3_median': c3_median
            })

        de_df = pd.DataFrame(de_results)

        # Apply Benjamini-Hochberg FDR correction
        # This adds a new column 'fdr'
        reject, pvals_corrected, _, _ = multipletests(de_df['p_value'], method='fdr_bh')
        de_df['fdr'] = pvals_corrected
        de_df.to_csv(DE_INFO_FILE, index=False, sep="\t")
    else:
        de_df = pd.read_csv(DE_INFO_FILE, sep='\t')
    # Filter: FDR < 0.05
    significant_genes_df = de_df[de_df['fdr'] < 0.05]
    de_genes_list = significant_genes_df['geneID'].unique()

    sig_with_coords = pd.merge(gene_coords_df[['chrom', 'start', 'end', 'geneID']],
                               significant_genes_df,
                               on='geneID',
                               how='inner')

    # Save Results
    sig_with_coords.to_csv('/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/DE_genes_FDR_0.05.tsv', index=False, sep="\t")

    print(f"Total Genes: {len(de_df)}")
    print(f"Significant DE Genes (FDR < 0.05): {len(de_genes_list)}")

    # ---------------------------------------------------------
    # 3. Enrichment Analysis (Updated with FDR list)
    # ---------------------------------------------------------
    print("\nCalculating Enrichment...")

    all_genes = set(de_df['geneID'])
    de_genes_set = set(de_genes_list)

    # Union of genes near hyper OR hypo DMRs
    set_hyper = set(genes_near_hyper)
    set_hypo = set(genes_near_hypo)
    near_dmr_set = set_hyper.union(set_hypo)

    # Contingency Table
    #              Near DMR    Not Near DMR
    # DE (FDR<.05) a           b
    # Not DE       c           d

    a = len(de_genes_set.intersection(near_dmr_set))
    b = len(de_genes_set - near_dmr_set)
    c = len((all_genes - de_genes_set).intersection(near_dmr_set))
    d = len((all_genes - de_genes_set) - near_dmr_set)

    if (a + b + c + d) > 0:
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
        print(f"Total DE Genes: {len(de_genes_set)}")
        print(f"DE Genes Near DMRs: {a}")
        print(f"DE Genes Not Near DMRs: {b}")
        print(f"Non-DE Genes Near DMRs: {c}")
        print(f"Non-DE Genes Not Near DMRs: {d}")
        print(f"Enrichment Odds Ratio: {oddsratio:.4f}")
        print(f"Fisher Exact p-value: {pvalue:.4e}")
    else:
        print("Insufficient data for enrichment analysis.")

    run_matched_enrichment(de_df, de_genes_set, near_dmr_set, k=10)
    # ---------------------------------------------------------
    # 4. Classification: Enhancers vs Silencers
    # ---------------------------------------------------------
    print("\nClassifying Enhancers and Silencers...")
    enhancer_analysis(df_dmr, genes_bed, significant_genes_df)



def plot_meth_vs_expression(df, output_path):
    """
    Generates scatter plot using the provided styling template.
    """
    print("Generating plot...")
    plt.figure(figsize=(10, 8))

    # Create the scatter plot
    # We map the calculated metrics to x and y.
    # We map the regulatory type to hue.
    sns.scatterplot(
        data=df,
        x='meth_diff_B6_minus_C3',
        y='log2FC_expr_B6_vs_C3',
        hue='Regulatory Type',
        palette={
            'Enhancer (Negative Correlation)': 'blue',
            'Silencer (Positive Correlation)': 'red'
        },
        alpha=0.7,
        s=80,
        legend=False  # Changed to True so we can see which color is which
    )

    # Add a horizontal line at y=0 (Equal expression)
    plt.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

    # Add a vertical line at x=0 (Equal methylation)
    plt.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # Styling aesthetic (from template)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add labels and title
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=18)

    # Updated axes labels to reflect the specific B6 vs C3 comparison
    plt.xlabel('Methylation Difference (B6 - C3)', fontsize=22)
    plt.ylabel('Expression Log2 Fold Change (B6 / C3)', fontsize=22)
    plt.title('DMR Methylation Diff vs Gene Expression FC', fontsize=24)

    # Add grid
    plt.grid(True, alpha=0.3)

    # Customize Legend options
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=12, title_fontsize=14)

    plt.tight_layout()
    print(f"Saving plot to {output_path}")
    # plt.savefig(output_path, dpi=300)
    plt.show() # Uncomment if running in an interactive environment like Jupyter


def load_and_prepare_data():
    """
    Loads enhancer and silencer results, combines them, and calculates
    plotting metrics (methylation difference vs expression fold change).
    """
    print("Loading data...")
    if not os.path.exists(ENHANCER_FILE) or not os.path.exists(SILENCER_FILE):
        print(f"Error: Input files not found in {BASE_DIR}. Please run previous analysis step first.")
        return None

    df_enh = pd.read_csv(ENHANCER_FILE, sep='\t')
    df_sil = pd.read_csv(SILENCER_FILE, sep='\t')

    if len(df_enh) == 0 and len(df_sil) == 0:
        print("No candidates found in input files. Cannot plot.")
        return None

    # Add labels for coloring
    df_enh['Regulatory Type'] = 'Enhancer (Negative Correlation)'
    df_sil['Regulatory Type'] = 'Silencer (Positive Correlation)'

    # Combine into one dataframe
    df = pd.concat([df_enh, df_sil], ignore_index=True)

    # --- Metrics Calculation ---

    # 1. X-axis: Methylation Difference (B6 - C3)
    # Interpretation: Positive values mean B6 is more methylated than C3.
    # Negative values mean B6 is less methylated (C3 is higher).
    df['meth_diff_B6_minus_C3'] = df['DMR_B6_meth'] - df['DMR_C3_meth']

    # 2. Y-axis: Expression Log2 Fold Change (B6 vs C3)
    # We add a small pseudo-count (+1) to avoid errors if median expression is 0.
    # Interpretation: Positive values mean B6 expression is higher.
    # Negative values mean C3 expression is higher.
    # Using log2 is standard for expression fold changes.
    df['log2FC_expr_B6_vs_C3'] = np.log2(
        (df['B6_median_expr'] + 1) / (df['C3_median_expr'] + 1)
    )

    print(f"Data prepared: {len(df)} total gene-DMR pairs.")
    return df


if __name__ == "__main__":
    load_and_analyze()
    df = load_and_prepare_data()
    plot_meth_vs_expression(df, "")