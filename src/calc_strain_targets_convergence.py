import pandas as pd
import scipy.stats as stats
import glob
import os
import math
from pybedtools import BedTool

# --- Configuration ---
GENE_COUNTS_FILE = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/GSE106208_counts_per_gene.tpm.txt.gz'
GENES_BED_MM10 = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/genes_mm10.bed'
DMR_FILE = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/B6_C3.DMRS.tsv'
TF_DIR = '/home/jon/Workspace/sd_asm_paper/review_comments/mouse_rna/mouse_tfs'
WINDOW_SIZE = 150000  # 150kb


def get_stats(set_a, set_b, universe_size):
    """
    Calculates overlap, expected overlap, fold change, and p-value between two sets.
    """
    len_a = len(set_a)
    len_b = len(set_b)
    overlap_count = len(set_a.intersection(set_b))

    # Expected Overlap (Mean of Hypergeometric)
    if universe_size > 0:
        expected = stats.hypergeom.mean(universe_size, len_a, len_b)
    else:
        expected = 0

    # Fold Change (Observed / Expected)
    if expected > 0:
        fold_change = overlap_count / expected
    else:
        fold_change = 0

    # Hypergeometric P-value
    p_val = stats.hypergeom.sf(overlap_count - 1, universe_size, len_a, len_b)

    return overlap_count, expected, fold_change, p_val


def load_and_analyze():
    print("Loading Global Data...")

    # 1. Load Genes (Universe)
    genes_bed = BedTool(GENES_BED_MM10)
    gene_coords_df = pd.read_csv(GENES_BED_MM10, sep='\t', header=None,
                                 names=['chrom', 'start', 'end', 'geneID', 'score'])
    N = gene_coords_df['geneID'].nunique()
    print(f"Universe Size (N): {N}")

    # 2. Load DMRs
    df_dmr = pd.read_csv(DMR_FILE, sep='\t')
    hyper_df = df_dmr[df_dmr['B6 status'] == 'hypermethyalted']
    hypo_df = df_dmr[df_dmr['B6 status'] == 'hypomethylated']

    # Convert to BedTool and Sort
    hyper_bed = BedTool.from_dataframe(hyper_df[['#chrom', 'start', 'end', 'region']]).sort()
    hypo_bed = BedTool.from_dataframe(hypo_df[['#chrom', 'start', 'end', 'region']]).sort()

    # ---------------------------------------------------------
    # Part 1: Define Global Sets (For the Enrichment Ratio Analysis)
    # ---------------------------------------------------------
    print("Identifying global gene sets...")

    global_genes_hyper = set(genes_bed.window(hyper_bed, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique())
    global_genes_hypo = set(genes_bed.window(hypo_bed, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique())
    ov_global, exp_global, fc_global, pval_global = get_stats(global_genes_hyper, global_genes_hypo, N)
    print("\n" + "=" * 60)
    print("GLOBAL ANALYSIS (All Genes, No TF Filtering)")
    print("=" * 60)
    print(f"Hyper Genes: {len(global_genes_hyper)}")
    print(f"Hypo Genes:  {len(global_genes_hypo)}")
    print(f"Overlap:     {ov_global}")
    print(f"Expected:    {exp_global:.2f}")
    print(f"Fold Change: {fc_global:.2f}")
    print(f"P-Value:     {pval_global:.5e}")
    print("=" * 60 + "\n")

    # ---------------------------------------------------------
    # Part 2: TF-Specific Analysis
    # ---------------------------------------------------------
    print("\n=== TF-Specific Analysis ===")

    tf_files = glob.glob(os.path.join(TF_DIR, "*.TF"))
    print(f"Found {len(tf_files)} TF files in {TF_DIR}")

    results = []

    for tf_file in tf_files:
        tf_name = os.path.basename(tf_file).replace('.TF', '')

        try:
            tf_bed = BedTool(tf_file).sort()
        except Exception as e:
            print(f"Skipping {tf_name}: Error reading file ({e})")
            continue

        # --- A. Filter DMRs by TF Intersection ---
        hyper_dmrs_tf = hyper_bed.intersect(tf_bed, u=True)
        hypo_dmrs_tf = hypo_bed.intersect(tf_bed, u=True)

        # Get Gene Sets for this TF
        if len(hyper_dmrs_tf) > 0:
            genes_hyper_tf = set(genes_bed.window(hyper_dmrs_tf, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique())
        else:
            genes_hyper_tf = set()

        if len(hypo_dmrs_tf) > 0:
            genes_hypo_tf = set(genes_bed.window(hypo_dmrs_tf, w=WINDOW_SIZE, u=True).to_dataframe()['name'].unique())
        else:
            genes_hypo_tf = set()

        # If this TF has no interactions at all, skip
        if len(genes_hyper_tf) == 0 and len(genes_hypo_tf) == 0:
            continue

        # --- B. Analysis 1: TF-Specific Overlap (Hyper_TF vs Hypo_TF) ---
        ov_tf, exp_tf, fc_tf, pval_tf = get_stats(genes_hyper_tf, genes_hypo_tf, N)

        # --- C. Analysis 2: Global Enrichment Ratio (TF_Set vs Global_Hyper / Global_Hypo) ---
        # "TF Set" = Genes near ANY DMR that has this TF binding site
        genes_near_tf = genes_hyper_tf.union(genes_hypo_tf)

        # Enrichment in Global Hyper
        obs_hyper, exp_hyper, fc_hyper, _ = get_stats(global_genes_hyper, genes_near_tf, N)

        # Enrichment in Global Hypo
        obs_hypo, exp_hypo, fc_hypo, _ = get_stats(global_genes_hypo, genes_near_tf, N)

        # Calculate Ratio
        if fc_hypo > 0:
            ratio = fc_hyper / fc_hypo
        else:
            ratio = 0

            # --- D. Store Results ---
        results.append({
            'TF': tf_name,
            # Data for Analysis 1
            'Hyper_Genes_TF': len(genes_hyper_tf),
            'Hypo_Genes_TF': len(genes_hypo_tf),
            'Overlap_TF': ov_tf,
            'Expected_TF': exp_tf,
            'FC_TF_Overlap': fc_tf,
            'P_Value_TF': pval_tf,

            # Data for Analysis 2
            'Obs_Hyper': obs_hyper,
            'Exp_Hyper': exp_hyper,
            'FC_Hyper': fc_hyper,
            'Obs_Hypo': obs_hypo,
            'Exp_Hypo': exp_hypo,
            'FC_Hypo': fc_hypo,
            'Ratio_Hyper_vs_Hypo': ratio
        })

    # ---------------------------------------------------------
    # Output Results
    # ---------------------------------------------------------
    if results:
        df_res = pd.DataFrame(results)

        # Display settings
        pd.set_option('display.max_rows', None)
        pd.set_option('display.width', 1000)

        # Formatters
        formatters = {
            'Expected_TF': '{:.2f}'.format,
            'FC_TF_Overlap': '{:.2f}'.format,
            'P_Value_TF': '{:.2e}'.format,
            'Exp_Hyper': '{:.2f}'.format,
            'FC_Hyper': '{:.2f}'.format,
            'Exp_Hypo': '{:.2f}'.format,
            'FC_Hypo': '{:.2f}'.format,
            'Ratio_Hyper_vs_Hypo': '{:.3f}'.format
        }

        # Columns to display
        cols = ['TF',
                'FC_TF_Overlap', 'Overlap_TF', 'Expected_TF', 'P_Value_TF',  # Analysis 1
                'Ratio_Hyper_vs_Hypo',
                'Obs_Hyper', 'Exp_Hyper',  # New: Hyper Enrichment details
                'Obs_Hypo', 'Exp_Hypo'  # New: Hypo Enrichment details
                ]

        # --- Print 1: Sorted by TF-Specific Overlap Fold Change ---
        print("\n" + "=" * 90)
        print("TABLE 1: Sorted by TF-Specific Overlap (FC_TF_Overlap)")
        print("Do the TF-associated Hyper and Hypo genes overlap more than expected?")
        print("=" * 90)
        df_sorted_overlap = df_res.sort_values(by='FC_TF_Overlap', ascending=False)
        print(df_sorted_overlap[cols].to_string(index=False, formatters=formatters))

        # --- Print 2: Sorted by Hyper/Hypo Ratio ---
        print("\n" + "=" * 90)
        print("TABLE 2: Sorted by Hyper vs Hypo Bias (Ratio)")
        print("Is the TF more strongly associated with Hyper genes or Hypo genes?")
        print("=" * 90)
        df_sorted_ratio = df_res.sort_values(by='Ratio_Hyper_vs_Hypo', ascending=False)
        print(df_sorted_ratio[cols].to_string(index=False, formatters=formatters))

    else:
        print("No results found.")


if __name__ == "__main__":
    load_and_analyze()