import pandas as pd
import os


if __name__ == '__main__':
    is_first = True
    count = 0
    input_dir= "sd_asm_analysis/homog/homog_aligned/all_snps/"
    res_df = pd.read_csv(os.path.join(input_dir, "big_hetero_homo_snps.tsv"), sep="\t")
    origi_cols = list(res_df.columns)
    res_df["num_hetero"] = (res_df[list(origi_cols[2:])] == 1).sum(axis=1)
    res_df["num_homo_a"] = (res_df[list(origi_cols[2:])] == 2).sum(axis=1)
    res_df["num_homo_b"] = (res_df[list(origi_cols[2:])] == 3).sum(axis=1)
    snps_with_many_hetero_low_homo = res_df[(res_df["num_hetero"] > 150) | (res_df["num_homo_a"] + res_df["num_homo_b"]  < 5)]
    snps_with_many_hetero_low_homo[["#chrom", "start", "num_hetero", "num_homo_a", "num_homo_b"]].to_csv(os.path.join(input_dir, "snps_with_many_hetero.tsv"), sep="\t", index=False)
    snps_with_homo_hetero = res_df[(res_df["num_hetero"] > 5) & (res_df["num_homo_a"] > 5) & (res_df["num_homo_b"] > 5)]
    snps_with_homo_hetero[["#chrom", "start", "num_hetero", "num_homo_a", "num_homo_b"]].to_csv(os.path.join(input_dir, "snps_with_homo_hetero.tsv"), sep="\t", index=False)
    x = 0
