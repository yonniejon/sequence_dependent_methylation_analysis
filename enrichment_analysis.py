#!/usr/bin/python3 -u
import argparse
from multiprocessing import Pool
from scipy.stats import norm
import pandas as pd
import uuid
from pathlib import Path
import os
import numpy as np

from asm_utils import run_command
from permute_blocks import permute


def get_intersection_command(count_blocks_in_this, intersecting_this, use_bp, use_bg, bed_graph_operation, count_anno):
    block_in_name, anno_in_name = ("b", "a") if count_anno else ("a", "b")

    if use_bg:
        if bed_graph_operation == "sum":
            awk_part = "awk 'BEGIN{sum=0}{sum+=$4}END{print sum/NR}'"
        elif bed_graph_operation == "mean":
            awk_part = "awk 'BEGIN{sum=0}{region_len=$3-$2 + 1; sum+=($4/(region_len))}END{print sum/NR}'"
        cmd = f"bedtools intersect -wa -{anno_in_name} {intersecting_this} -{block_in_name} {count_blocks_in_this}" \
              f" | sort -u | {awk_part}"
        # cmd = f"bedtools intersect -wa -a {count_blocks_in_this} -b {intersecting_this} | sort -u | {awk_part}"
        return cmd
    else:
        use_bp_suffix = "| sort -u | awk 'BEGIN{sum=0}{sum+=$NF}END{print sum}'"
        non_use_bp_suffix = "| sort -k1,1 -k2,2n | bedtools merge | wc -l "
        return f"bedtools intersect -{'wo' if use_bp else 'wa'} -{block_in_name} {count_blocks_in_this}" \
               f" -{anno_in_name} {intersecting_this} {use_bp_suffix if use_bp else non_use_bp_suffix}"
        # return f"bedtools intersect -{'wo' if use_bp else 'wa'} -b {intersecting_this} -a {count_blocks_in_this}
        # {use_bp_suffix if use_bp else non_use_bp_suffix}"


def get_intersection(block_file, annotation_file, use_bp, use_bg, bed_graph_operation, count_anno):
    intersection = run_command(get_intersection_command(block_file, annotation_file, use_bp, use_bg,
                                                        bed_graph_operation, count_anno)).strip()
    # try:
    #     int(intersection)
    # except:
    #     aba = 0
    return float(intersection)


def permute_save_and_delete(block_file, annotation_file, chrom_size_file, tmp_dir, use_bp, use_bg, bed_graph_operation,
                            count_anno):
    tmp_file = f"permuted_annotation_file_{str(uuid.uuid4())}.bed"
    tmp_file = os.path.join(tmp_dir, tmp_file)
    try:
        cur_seed = int.from_bytes(os.urandom(4), 'big')
        permuted_df = permute(block_file, chrom_size_file, seed=cur_seed)
        permuted_df.to_csv(tmp_file, sep='\t', index=False, header=False)
        intersection = get_intersection(tmp_file, annotation_file, use_bp, use_bg, bed_graph_operation, count_anno)
        return intersection
    finally:
        os.remove(tmp_file)


def get_random_blocks_from_background(background_blocks, n_rows, seed=42):
    block_df = pd.read_csv(background_blocks, sep='\t', usecols=[0, 1, 2], names=["chrom", "start", "end"])
    return block_df.sample(n=n_rows, random_state=seed)


def sample_and_get_intersection(background_file, annotation_file, n_rows, tmp_dir, use_bp, use_bg, bedgraph_op,
                                count_anno):
    tmp_file = f"permuted_annotation_file_{str(uuid.uuid4())}.bed"
    tmp_file = os.path.join(tmp_dir, tmp_file)
    try:
        cur_seed = int.from_bytes(os.urandom(4), 'big')
        sampled_df = get_random_blocks_from_background(background_file, n_rows, seed=cur_seed)
        sampled_df.to_csv(tmp_file, sep='\t', index=False, header=False)
        intersection = get_intersection(tmp_file, annotation_file, use_bp, use_bg, bedgraph_op, count_anno)
        return intersection
    finally:
        os.remove(tmp_file)


def main(args):
    tmp_dir = args.tmp_dir#"enrichment_analysis_temp_files"
    input_block_file = args.block_file #"/cs/zbio/jrosensk/imprinting_files/all_snps_dir/bimodal_with_poo_asm.bed"
    annotation_file = args.annotation_file# "/cs/cbio/jon/bio_annot/ctcfdb.merged.bed.gz"
    chrom_size_file = args.chrom_size_file#"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/chrome.size"
    num_iterations = int(args.num_iterations)
    num_threads = int(args.nr_threads)#1
    original_intersection = get_intersection(input_block_file, annotation_file, args.use_bp, args.use_bedgraph_vals,
                                             args.bed_graph_operation, args.count_anno)
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    if args.background_file is None:

        params = [(input_block_file, annotation_file, chrom_size_file, tmp_dir, args.use_bp, args.use_bedgraph_vals,
                   args.bed_graph_operation, args.count_anno) for _ in range(num_iterations)]
        with Pool(num_threads) as p:
            res = p.starmap(permute_save_and_delete, params)
            p.close()
            p.join()
    else:
        block_df = pd.read_csv(input_block_file, sep='\t', usecols=[0, 1, 2], names=["chrom", "start", "end"])
        num_blocks = block_df.shape[0]
        del block_df

        params = [(args.background_file, annotation_file, num_blocks, tmp_dir, args.use_bp, args.use_bedgraph_vals,
                   args.bed_graph_operation, args.count_anno) for _ in
                  range(num_iterations)]
        with Pool(num_threads) as p:
            res = p.starmap(sample_and_get_intersection, params)
            p.close()
            p.join()


    sample_mean = sum(res) / len(res)
    std_dev = np.std(res)
    estimated_p_val = norm.sf(original_intersection, loc=sample_mean, scale=std_dev)
    if original_intersection < sample_mean:
        estimated_p_val = 1 - estimated_p_val

    estimated_p_val = min(estimated_p_val * 2, 1)
    # count_above = sum([int(cur_intersection > original_intersection) for cur_intersection in res])
    # for i in range(num_iterations):
    #     cur_intersection = permute_save_and_delete(input_block_file, annotation_file, chrom_size_file)
    #     count_above += int(cur_intersection > original_intersection)
    # empirical_p = count_above / num_iterations
    print(f"Original Intersection: {original_intersection}")
    print(f"Permuted Intersections: {','.join([str(r) for r in res])}")
    print(f"Empirical p-value for enrichment: {estimated_p_val}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('block_file', help='Bed file')
    parser.add_argument('annotation_file', help='Bed file')
    parser.add_argument('chrom_size_file', help='File of sizes of chromosomes')
    parser.add_argument('--tmp_dir', default="enrichment_analysis_temp_files")
    parser.add_argument('--nr_threads', default=1)
    parser.add_argument('--num_iterations', default=10)
    parser.add_argument('--use_bp', help="count the number of base pairs intersecting, rather the number of"
                                         " intersecting annotation regions", action="store_true")
    parser.add_argument('--use_bedgraph_vals', action="store_true")
    parser.add_argument('--bed_graph_operation', default="sum", help="Options are sum or mean")
    parser.add_argument('--background_file', default=None)
    parser.add_argument('--count_anno', action="store_true")
    args = parser.parse_args()
    main(args)
