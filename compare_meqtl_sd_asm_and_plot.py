import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import os


def load_tissue_data(data_dir, bg_file, ov_file):
    """
    Loads background and overlap files, separates them into
    'Inside SD-ASM' and 'Outside SD-ASM' groups.
    """
    bg_file = os.path.join(data_dir, bg_file)
    ov_file = os.path.join(data_dir, ov_file)
    print(f"Processing: {os.path.basename(bg_file)}")

    if not os.path.exists(bg_file) or not os.path.exists(ov_file):
        print(f"  ERROR: File not found ({bg_file} or {ov_file})")
        return None, None

    # Load Data (Cols: 0=Chr, 1=Start, 3=Beta)
    # Assumes input BED format: Chr, Start, End, Beta, SNP_ID
    df_bg = pd.read_csv(bg_file, sep="\t", header=None, usecols=[0, 1, 3])
    df_bg.columns = ['chr', 'start', 'beta']

    df_ov = pd.read_csv(ov_file, sep="\t", header=None, usecols=[0, 1, 3])
    df_ov.columns = ['chr', 'start', 'beta']

    # Create Unique Keys (Chr:Start) to identify CpGs
    df_bg['key'] = df_bg['chr'].astype(str) + ":" + df_bg['start'].astype(str)
    df_ov['key'] = df_ov['chr'].astype(str) + ":" + df_ov['start'].astype(str)

    # Define Sets
    overlap_keys = set(df_ov['key'])

    # Group 1: Inside SD-ASM
    # (Drop duplicates in case one CpG overlaps multiple ASM regions)
    betas_inside = df_ov.drop_duplicates(subset='key')['beta']

    # Group 2: Outside SD-ASM (Background MINUS Overlap)
    df_outside = df_bg[~df_bg['key'].isin(overlap_keys)]
    betas_outside = df_outside['beta']

    return betas_inside, betas_outside


def run_stats(betas_in, betas_out, tissue_name):
    """Calculates Medians and Mann-Whitney U P-value"""
    print(f"\n--- Statistics for {tissue_name} ---")
    n_in = len(betas_in)
    n_out = len(betas_out)

    med_in = betas_in.median()
    med_out = betas_out.median()
    mean_in = betas_in.mean()
    mean_out = betas_out.mean()

    print(f"  Count (Inside):  {n_in}")
    print(f"  Count (Outside): {n_out}")
    print(f"  Median Beta (Inside):  {med_in:.4f}")
    print(f"  Median Beta (Outside): {med_out:.4f}")
    print(f"  Mean Beta (Inside):  {mean_in:.4f}")
    print(f"  Mean Beta (Outside): {mean_out:.4f}")

    # Test if Inside > Outside
    stat, p_val = stats.mannwhitneyu(betas_in, betas_out, alternative='greater')
    print(f"  P-value: {p_val:.5e}")

    return p_val, med_in, med_out


def main():
    # ================= FILES CONFIGURATION =================
    # Update these filenames if they differ on your system
    data_dir = "/cs/cbio/jon/sd_asm_hg38/reviewer_comments/jaffe/jaffe_meqtls"
    files = {
        'DLPFC': {
            'bg': 'dlpfc_mCpGs.150_dist.bed',
            'ov': 'dlpfc_mcpgs.neuron_sd_asm.intersected.bed'
        },
        'Hippocampus': {
            'bg': 'hippo_mCpGs.150_dist.bed',
            'ov': 'hippo_mcpgs.neuron_sd_asm.intersected.bed'
        }
    }
    OUTPUT_PLOT = "meqtl_sdasm_comparison_combined.png"
    # =======================================================

    # Prepare plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    # Loop through tissues
    for i, (tissue, paths) in enumerate(files.items()):
        betas_in, betas_out = load_tissue_data(data_dir, paths['bg'], paths['ov'])

        if betas_in is None:
            continue

        # Run Stats
        p_val, med_in, med_out = run_stats(betas_in, betas_out, tissue)

        # Plotting
        ax = axes[i]
        data_to_plot = [betas_out, betas_in]

        # Create boxplot
        bp = ax.boxplot(data_to_plot, showfliers=False, patch_artist=True, widths=0.6)

        # Colors
        colors = ['#d9d9d9', '#6baed6']  # Grey for Outside, Blue for Inside
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)

        # Labels and Title
        ax.set_xticklabels(['Outside\nSD-ASM', 'Inside\nSD-ASM'], fontsize=10)
        ax.set_title(f"{tissue}\n(MWU p={p_val:.1e})", fontsize=12, fontweight='bold')
        ax.grid(axis='y', linestyle='--', alpha=0.5)

        if i == 0:
            ax.set_ylabel('Max Absolute meQTL Beta', fontsize=11)

    plt.suptitle("Impact of SD-ASM on Local meQTL Effect Size", fontsize=14)
    plt.tight_layout()
    plt.show()
    # plt.savefig(OUTPUT_PLOT, dpi=300)
    # print(f"\nPlot saved to {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()