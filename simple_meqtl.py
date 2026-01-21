import os
from multiprocessing import Pool
import pandas as pd
import numpy as np
import optparse
import scipy.stats as stats
from fdr_corection_allele_homog import choose_blocks_by_fdr_bh
from asm_utils import format_float
from handle_snps import classify_as_imprinting
import gc

# FILENAME_SUFFIX = 'allele.homog.pval.gz'
FILENAME_SUFFIX = 'allele.homog.gz'


def classify_as_sd(num_times_asm, num_times_asm_let_1_u):
    return num_times_asm >= 1 and (0 == num_times_asm_let_1_u or num_times_asm_let_1_u == num_times_asm)


def iterate_over_files_with_pandas(file_name, base_name, with_p_val=True, filter_non_asm=False, num_rows=None, female_x=False):
    if with_p_val:
        expected_cols = ["chrom", "start", "let1", "let2", "u1", "x1", "m1", "u2", "x2", "m2", "p_val", "let1_is_u"]
        final_expected_cols = expected_cols[:-2]
    else:
        expected_cols = ["chrom", "start", "let1", "let2", "u1", "x1", "m1", "u2", "x2", "m2"]

    zipped = list(zip(expected_cols, range(0, len(expected_cols))))
    zipped = zipped[0:2] + zipped[-2:]
    if num_rows is None:
        cur_homog = pd.read_csv(file_name, sep="\t", names=expected_cols)
    else:
        cur_homog = pd.read_csv(file_name, sep="\t", names=expected_cols, nrows=num_rows)
    if female_x:
        cur_homog = cur_homog[cur_homog["chrom"] == "chrX"]

    if filter_non_asm:
        cur_homog = cur_homog[cur_homog["p_val"] >= 0]
    if with_p_val:
        cur_homog = cur_homog[final_expected_cols]
    cur_homog["name"] = base_name
    return cur_homog


def meqtl_test_agg(row, num_samples_included, psuedo_count=0, high_threshold=0.6, low_threshold=0.4,
                   min_read_count_per_allele=5, require_u_higher_dif_alleles=True,
                   require_hetero_and_homo=False, hetero_homo_thresh=2, female_x=False):
    counts_array = np.zeros([2, 2])
    counts_array[0, 0] = row['u1']
    counts_array[0, 1] = row['m1']
    counts_array[1, 0] = row['u2']
    counts_array[1, 1] = row['m2']
    counts_array = counts_array + psuedo_count

    total_1 = row['u1'] + row['m1'] + row['x1']
    total_2 = row['u2'] + row['m2'] + row['x2']
    prop_1 = row['u1'] / total_1 if total_1 > 0 else 0
    prop_2 = row['u2'] / total_2 if total_2 > 0 else 0
    avg_read_count_1 = total_1 / num_samples_included
    avg_read_count_2 = total_2 / num_samples_included

    eps = 0.00001
    res_letter = 'NA'
    p = 1
    if (avg_read_count_1 > min_read_count_per_allele and avg_read_count_2 > min_read_count_per_allele):
        u_higher_in_1_low_2 = row['u1'] / (row['m1'] + eps) > 1 and row['u2'] / (row['m2'] + eps) < 1
        u_higher_in_2_low_1 = row['u2'] / (row['m2'] + eps) > 1 and row['u1'] / (row['m1'] + eps) < 1
        if (u_higher_in_2_low_1 or u_higher_in_1_low_2 or female_x):
            oddsr, p = stats.fisher_exact(counts_array, alternative="two-sided")

    return p


def iterate_over_groups(in_file=None):
    if in_file is None:
        group_file = "groups207.csv"
    else:
        group_file = in_file
    gdf = pd.read_csv(group_file)
    gdf = gdf[gdf['include'] == True]
    # gdf = gdf[gdf['source'] == 'Dor']
    mgroups = list(sorted(gdf['group'].unique()))
    groups = [g for g in mgroups if ':' not in g]
    group_dict = {}
    for group in groups:
        samples = sorted(gdf[gdf['group'] == group]['name'])
        group_files = []
        for sample in samples:
            cur_file_name = sample + "." + FILENAME_SUFFIX
            group_files.append(cur_file_name)
        group_dict[group] = set(group_files)
    return group_dict


def get_snp_dict_tissue_specific(tissue_name, input_dir, num_threads, significance_threshold, require_u_higher_dif_alleles_in_test,
                                 outdir, female_x):
    group_dict = iterate_over_groups(None)
    require_hetero_and_homo_in_test = False
    num_rows = None
    filter_non_asm = False
    global p_val_list, let1_is_u_list, snp_list, snp_dict, snp, cur_dict, sample_name_list
    is_first = True
    file_list = []
    group_p_vals = []
    count = 0
    for group_name, samples in list(group_dict.items()):
        out_file = os.path.join(outdir, f"{group_name}.sd_asm_chart.tsv")
        if female_x:
            out_file = os.path.join(outdir, f"{group_name}.chrX.sd_asm_chart.tsv")
        if os.path.isfile(out_file):
            print(f"{out_file} exists")
            continue
        print(out_file)

        samples = list(samples)

        file_list = []
        for s in samples:
            file_name = os.path.join(input_dir, s) # + "." + FILENAME_SUFFIX
            base_name = s.split(".")[0]
            if os.path.isfile(file_name):
                count += 1
                file_list.append((os.path.join(input_dir, s), base_name, False, filter_non_asm, num_rows, female_x))
        print(group_name)
        with Pool(num_threads) as p:
            all_dfs = p.starmap(iterate_over_files_with_pandas, file_list)
            p.close()
            p.join()
        if len(all_dfs) < 2:
            continue
        res_df = pd.concat(all_dfs)

        snp_list = res_df["chrom"] + ":" + res_df["start"].apply(str)

        grouped_res_df = res_df.groupby(['chrom', 'start']).agg(
            {"u1": 'sum', "x1": 'sum', "m1": 'sum', "u2": 'sum', "x2": 'sum', "m2": 'sum'})

        del res_df
        grouped_res_df['meqtl_p_val'] = grouped_res_df.apply(lambda row: meqtl_test_agg(row, len(file_list),
                                                                                        require_u_higher_dif_alleles=require_u_higher_dif_alleles_in_test,
                                                                                        require_hetero_and_homo=require_hetero_and_homo_in_test, female_x=female_x),
                                                             axis=1)
        correct_p_vals, all_snps_and_sample_names, _ = choose_blocks_by_fdr_bh(list(grouped_res_df['meqtl_p_val']),
                                                                               list(grouped_res_df['meqtl_p_val']),
                                                                               alpha=significance_threshold,
                                                                               return_only_selected=False)
        grouped_res_df['corrected_meqtl_p_val'] = correct_p_vals
        psuedo_count = 1
        high_thresh = 0.7
        grouped_res_df = grouped_res_df.reset_index()
        grouped_res_df['u_prop_1'] = (grouped_res_df['u1'] + psuedo_count) / (
                    grouped_res_df['u1'] + grouped_res_df['m1'] + grouped_res_df['x1'] + psuedo_count)
        grouped_res_df['u_prop_2'] = (grouped_res_df['u2'] + psuedo_count) / (
                grouped_res_df['u2'] + grouped_res_df['m2'] + grouped_res_df['x2'] + psuedo_count)
        grouped_res_df['m_prop_1'] = (grouped_res_df['m1'] + psuedo_count) / (
                grouped_res_df['u1'] + grouped_res_df['m1'] + grouped_res_df['x1'] + psuedo_count)
        grouped_res_df['m_prop_2'] = (grouped_res_df['m2'] + psuedo_count) / (
                grouped_res_df['u2'] + grouped_res_df['m2'] + grouped_res_df['x2'] + psuedo_count)
        grouped_res_df['is_u'] = (
                    (grouped_res_df['u_prop_1'] > high_thresh) & (grouped_res_df['u_prop_2'] > high_thresh)
                    & (grouped_res_df['corrected_meqtl_p_val'] >= significance_threshold))
        grouped_res_df['is_m'] = (grouped_res_df['m_prop_1'] > high_thresh) & (
                    grouped_res_df['m_prop_2'] > high_thresh) & (
                                         grouped_res_df['corrected_meqtl_p_val'] >= significance_threshold)
        grouped_res_df['is_b'] = (~grouped_res_df['is_u']) & (~grouped_res_df['is_m']) & (
                    grouped_res_df['corrected_meqtl_p_val'] < significance_threshold)
        grouped_res_df['is_x'] = (~grouped_res_df['is_u']) & (~grouped_res_df['is_m']) & (
                grouped_res_df['corrected_meqtl_p_val'] >= significance_threshold)
        grouped_res_df["status"] = "NA"
        grouped_res_df.loc[grouped_res_df["is_u"], 'status'] = "U"
        grouped_res_df.loc[grouped_res_df["is_m"], 'status'] = "M"
        grouped_res_df.loc[grouped_res_df["is_b"], 'status'] = "B"
        grouped_res_df.loc[grouped_res_df["is_x"], 'status'] = "X"
        grouped_res_df.loc[((grouped_res_df['u1'] + grouped_res_df['m1'] + grouped_res_df['x1'] < 10) |
                            (grouped_res_df['u2'] + grouped_res_df['m2'] + grouped_res_df['x2'] < 10)), 'status'] = "NA"

        grouped_res_df.loc[
            ((grouped_res_df["is_b"]) & (grouped_res_df['corrected_meqtl_p_val'] > significance_threshold)), 'status'] = "X"

        grouped_res_df = grouped_res_df[
            ["chrom", "start", "u1", "x1", "m1", "u2", "x2", "m2", "meqtl_p_val", "corrected_meqtl_p_val", "status"]]
        grouped_res_df = grouped_res_df.rename(columns={'chrom': '#chrom'})

        grouped_res_df.to_csv(out_file, sep="\t", index=False)
        del grouped_res_df
        gc.collect()


def add_all_files(input_dir, female_x):
    is_first = True
    count = 0
    group_cols = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".tsv.gz") and not "full_sd_asm" in filename:#.gz
            if female_x and not "chrX" in filename:
                continue
            base_name = filename.split(".")[0]
            full_file = os.path.join(input_dir, filename)
            df = pd.read_csv(full_file, sep="\t")
            if is_first:
                res_df = df[["#chrom", "start"]]
                is_first = False
            res_df[base_name] = df["status"]
            group_cols.append(base_name)
            count += 1
            print(f"{base_name}, {count}")
    if female_x:
        res_df.fillna("N").to_csv(os.path.join(input_dir, "full_sd_asm_all_cell_types.chrX.tsv"),
            sep="\t", index=False)  # full_sd_asm_all_cell_types
    else:
        res_df.fillna("N").to_csv(os.path.join(input_dir, "full_sd_asm_all_cell_types.tsv"),
                                  sep="\t", index=False)  # full_sd_asm_all_cell_types


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--tissue_name',
                      default="")
    parser.add_option('--num_threads',
                      default=16)
    parser.add_option('--female_chrX', action="store_true")
    options, arguments = parser.parse_args()
    tissue_name = options.tissue_name

    input_dir = "sd_asm_analysis/homog/homog_aligned/all_snps/"
    num_threads = int(options.num_threads)
    # tissue_specific_flag = False
    outdir = "sd_asm_analysis/homog/homog_aligned/all_snps/meqtl/simple_sd_asm/"
    get_snp_dict_tissue_specific(tissue_name, input_dir, num_threads, 0.05, False, outdir, options.female_chrX)
    add_all_files("sd_asm_analysis/homog/homog_aligned/all_snps/meqtl/simple_sd_asm",
                  options.female_chrX)
    # a = 0
