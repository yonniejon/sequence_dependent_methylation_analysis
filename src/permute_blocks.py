#!/usr/bin/python3 -u
import argparse
import pandas as pd
import numpy as np


def permute(original_block_file, chrom_size_file, seed=3797736516):
    np.random.seed(seed)
    block_df = pd.read_csv(original_block_file, sep='\t', usecols=[0, 1, 2], names=["chrom", "start", "end"])
    chrom_df = pd.read_csv(chrom_size_file, sep='\t', usecols=[0, 1], names=["chrom", "end_index"])
    block_df.sort_values(by=['chrom', 'start'], ignore_index=True, inplace=True)
    all_chroms = set(block_df.chrom)
    # all_chroms = ["chr16"]
    df_res_list = []
    for chrom in all_chroms:
        cur_end = list(chrom_df[chrom_df.chrom == chrom].end_index)
        if len(cur_end) == 0:
            print(f"Chromsome {chrom} is not in chrom_size_file. Skipping")
            continue
        cur_end = cur_end[0]
        cur_df = block_df[block_df.chrom == chrom]
        block_lengths = np.array(cur_df.end - cur_df.start)
        background_starts = np.array([0] + list(cur_df.end))
        background_ends = np.array(list(cur_df.start) + [cur_end])
        background_lengths = background_ends - background_starts
        if background_lengths[-1] < 0:
            print(f"Warning: {chrom} contains rows after chromosome end. Ignoring these rows.")
            cur_df = cur_df[(cur_df["start"] < cur_end) & (cur_df["end"] < cur_end)]
            block_lengths = np.array(cur_df.end - cur_df.start)
            background_starts = np.array([0] + list(cur_df.end))
            background_ends = np.array(list(cur_df.start) + [cur_end])
            background_lengths = background_ends - background_starts
        if sum(background_lengths < 0) > 0:
            print(f"Error chromosome {chrom}: Make sure input blocks file {original_block_file} does not have overlapping blocks and that chrome sizes are from the reference")
        np.random.shuffle(background_lengths)
        np.random.shuffle(block_lengths)
        background_status = np.zeros(background_lengths.shape)
        block_status = np.ones(block_lengths.shape)
        all_lengths = np.empty((background_lengths.size + block_lengths.size,), dtype=block_lengths.dtype)
        all_status = np.empty((background_status.size + block_status.size,), dtype=block_status.dtype)
        all_lengths[0::2] = background_lengths
        all_lengths[1::2] = block_lengths
        all_status[0::2] = background_status
        all_status[1::2] = block_status
        new_ends = np.cumsum(all_lengths) + 1
        new_starts = new_ends - all_lengths
        chrom_list = np.array([chrom] * len(all_lengths))
        new_df = pd.DataFrame({"chrom": chrom_list, "start": new_starts, "end": new_ends, "status": all_status})
        new_df = new_df[new_df.status == 1]
        new_df = new_df.drop('status', axis=1)
        df_res_list.append(new_df)
    final_df = pd.concat(df_res_list)
    return final_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('block_file', help='Bed file')
    parser.add_argument('chrom_size_file', help='File of sizes of chromosomes')
    parser.add_argument('out_file', help='The name of the output file.')
    args = parser.parse_args()

    final_df = permute(args.block_file, args.chrom_size_file)

    final_df.to_csv(args.out_file, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()