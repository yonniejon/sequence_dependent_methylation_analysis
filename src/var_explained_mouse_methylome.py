import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os


# ---------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------

def parse_metadata(sample_names):
    """
    Parses sample names to determine Strain and Tissue Type based on user rules.
    """
    meta = []
    for sample in sample_names:
        strain = "Unknown"
        tissue_type = "Liver"  # Default assumption

        # Rule 1: Determine Tissue
        if "mouse_" in sample:
            tissue_type = "Other Tissue"
        else:
            tissue_type = "Liver"

        # Rule 2: Determine Strain
        if "B6C3F1" in sample or "C3B6F1" in sample:
            strain = "Cross"
        elif "_C3_" in sample:
            strain = "C3"
        elif "_B6_" in sample:
            strain = "B6"
        elif tissue_type == "Other Tissue":
            strain = "B6"

        meta.append({
            "SampleID": sample,
            "Strain": strain,
            "Tissue_Type": tissue_type
        })

    return pd.DataFrame(meta)


def get_symmetric_limits(coords):
    """
    Calculates a common min and max range for both axes to ensure
    the plot is square and uses the same scale.
    """
    val_min = np.min(coords)
    val_max = np.max(coords)

    # Add 10% padding for better spacing
    padding = (val_max - val_min) * 0.1
    limit_min = val_min - padding
    limit_max = val_max + padding

    return limit_min, limit_max


def setup_square_plot(title, x_label, y_label, limits):
    """
    Creates a figure with fixed axes dimensions to ensure all plots
    are exactly the same size regardless of tick label length.
    """
    # Create figure
    fig = plt.figure(figsize=(8, 8))

    # Manually define axes position [left, bottom, width, height] (0 to 1 scale)
    # This forces the box to be identical in every plot
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])

    # Set titles and labels (Labelsize 15 for axis title, 11 for ticks)
    ax.set_title(title, fontsize=15)
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel(y_label, fontsize=15)

    # Set tick label size to 11
    ax.tick_params(axis='both', which='major', labelsize=11)

    # Force Square Limits
    ax.set_xlim(limits)
    ax.set_ylim(limits)

    # Force Square Aspect Ratio
    ax.set_aspect('equal', adjustable='box')

    return fig, ax


# ---------------------------------------------------------
# Analysis Functions
# ---------------------------------------------------------

def run_pca_liver_only(df, meta_df):
    print("Running PCA on Liver samples...")
    liver_samples = meta_df[meta_df["Tissue_Type"] == "Liver"]["SampleID"]
    X_liver = df[liver_samples].T
    meta_liver = meta_df[meta_df["SampleID"].isin(liver_samples)]

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_liver)
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    l_min, l_max = get_symmetric_limits(coords)

    # Setup Plot
    fig, ax = setup_square_plot(
        title="PCA: Liver Samples Only",
        x_label=f"PC1 ({pc1_var:.1f}%)",
        y_label=f"PC2 ({pc2_var:.1f}%)",
        limits=(l_min, l_max)
    )

    sns.scatterplot(
        x=coords[:, 0], y=coords[:, 1],
        hue=meta_liver["Strain"],
        palette={"C3": "red", "B6": "blue", "Cross": "purple"},
        s=150, edgecolor='black', legend=False, ax=ax
    )

    output_file = "/home/jon/Workspace/sd_asm_paper/review_comments/var_explained/PCA_Liver_Only.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved {output_file}")
    plt.close()


def run_pca_all_samples(df, meta_df):
    print("Running PCA on All samples...")
    sample_cols = meta_df["SampleID"].values
    X_all = df[sample_cols].T

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_all)
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    meta_df = meta_df.copy()
    meta_df["PC1"] = coords[:, 0]
    meta_df["PC2"] = coords[:, 1]

    l_min, l_max = get_symmetric_limits(coords)

    # Setup Plot
    fig, ax = setup_square_plot(
        title="PCA: All Samples",
        x_label=f"PC1 ({pc1_var:.1f}%)",
        y_label=f"PC2 ({pc2_var:.1f}%)",
        limits=(l_min, l_max)
    )

    sns.scatterplot(
        data=meta_df, x="PC1", y="PC2",
        hue="Strain", style="Tissue_Type",
        palette={"C3": "red", "B6": "blue", "Cross": "purple"},
        markers={"Liver": "o", "Other Tissue": "X"},
        s=150, edgecolor='black', legend=False, ax=ax
    )

    output_file = "/home/jon/Workspace/sd_asm_paper/review_comments/var_explained/PCA_All_Samples.png"
    plt.savefig(output_file, dpi=300)
    print(f"Saved {output_file}")
    plt.close()


# ---------------------------------------------------------
# Main Execution
# ---------------------------------------------------------

def main():
    input_file = "/home/jon/Workspace/sd_asm_paper/review_comments/var_explained/avg_meth.no_NA.low_var_filtered.point_1_min_var.table"

    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    print(f"Loading {input_file}...")
    df = pd.read_csv(input_file, delim_whitespace=True)

    metadata_cols = ['chr', 'start', 'end', 'startCpG', 'endCpG']
    sample_cols = [c for c in df.columns if c not in metadata_cols]

    meta_df = parse_metadata(sample_cols)

    # Run Analysis
    run_pca_liver_only(df, meta_df)
    run_pca_all_samples(df, meta_df)

    print("All plots generated successfully.")


if __name__ == "__main__":
    main()