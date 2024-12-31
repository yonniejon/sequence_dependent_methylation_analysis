import argparse
import gzip as gz
import os
import warnings


def get_boundary(boundary_it):
    cur_start_boundary_line = boundary_it.readline()
    start_boundary_tokens = cur_start_boundary_line.split("\t")
    return start_boundary_tokens[0], int(start_boundary_tokens[1])


def get_chrom_num(chrom1):
    if chrom1[3:].isdigit():
        return int(chrom1[3:])
    else:
        return -1


def find_new_boundary(boundary_it, cand_index, cur_position, prev_position, prev_chrom, cur_chrom, block_chrom):
    cur_start_bound_distance = abs(cand_index - cur_position)
    prev_start_bound_distance = abs(cand_index - prev_position)
    prev_chrom_num = get_chrom_num(prev_chrom)
    block_chrom_num = get_chrom_num(block_chrom)
    prev_prev_position = None

    while cur_start_bound_distance < prev_start_bound_distance or ((prev_chrom_num > 0 and block_chrom_num > 0) and (prev_chrom_num < block_chrom_num)):
        prev_prev_position = prev_position
        prev_position = cur_position
        prev_chrom = cur_chrom
        prev_chrom_num = get_chrom_num(prev_chrom)
        cur_chrom, cur_position = get_boundary(boundary_it)
        # if chrom_of_boundary != cur_chrom:
        #     break
        cur_start_bound_distance = abs(cand_index - cur_position)
        prev_start_bound_distance = abs(cand_index - prev_position)
    return prev_position, cur_position, prev_chrom, cur_chrom, prev_prev_position


def block_boundary_iterator(block_file, boundary_file, out_file, force_delete):
    if os.path.isfile(out_file):
        if force_delete:
            warnings.warn(f"out_file {out_file} exists. Deleting!")
            os.remove(out_file)
        else:
            print(f"out_file {out_file} exists. Exiting. To delete previous file use -f option")
            return
    with gz.open(block_file, 'rt') as block_lines, open(boundary_file, 'rt') as start_boundaries, open(boundary_file,
                                                                                                       'rt'
                                                                                                       ) as end_boundaries:
        prev_chrom, prev_start_boundary = get_boundary(start_boundaries)
        cur_chrom, cur_start_boundary = get_boundary(start_boundaries)
        get_boundary(end_boundaries)
        prev_end_chrom, prev_end_boundary = get_boundary(end_boundaries)
        cur_end_chrom, cur_end_boundary = get_boundary(end_boundaries)
        prev_prev_start_boundary = None
        has_next = True
        to_print_list = []
        counter = 0
        while has_next:
            try:
                cur_block = block_lines.readline()
                if cur_block:
                    if cur_block.startswith("#"):
                        continue
                    counter += 1
                    block_tokens = cur_block.split("\t")
                    block_start = int(block_tokens[3])
                    block_end = int(block_tokens[4])
                    block_chrom = block_tokens[0]
                    while block_start <= prev_end_boundary and has_next:
                        cur_block = block_lines.readline()
                        if cur_block:
                            block_tokens = cur_block.split("\t")
                            block_start = int(block_tokens[3])
                            block_end = int(block_tokens[4])
                            block_chrom = block_tokens[0]
                        else:
                            has_next = False
                            break
                    if not has_next:
                        break
                    prev_start_boundary, cur_start_boundary, prev_chrom, cur_chrom, prev_prev_start_boundary = find_new_boundary(
                        start_boundaries,
                        block_start,
                        cur_start_boundary,
                        prev_start_boundary,
                        prev_chrom,
                        cur_chrom,
                        block_chrom)
                    prev_end_boundary, cur_end_boundary, prev_chrom, cur_chrom, _ = find_new_boundary(end_boundaries,
                                                                                                      block_end,
                                                                                                      cur_end_boundary,
                                                                                                      prev_end_boundary,
                                                                                                      prev_chrom,
                                                                                                      cur_chrom,
                                                                                                      block_chrom)
                    to_print_start = prev_start_boundary
                    to_print_end = prev_end_boundary
                    if prev_start_boundary == prev_end_boundary:
                        if block_start <= prev_start_boundary and block_end > prev_end_boundary:
                            continue
                        else:
                            if block_start < prev_start_boundary:
                                if prev_prev_start_boundary is not None:
                                    to_print_start = prev_prev_start_boundary
                                else:
                                    warnings.warn(f"Was not able to find appropriate boundaries for block {block_start} {block_end}")
                                    continue
                            else:
                                prev_end_boundary = cur_end_boundary
                                cur_end_chrom, cur_end_boundary = get_boundary(end_boundaries)
                                to_print_end = prev_end_boundary

                    to_print = "{}\t{}\n".format(to_print_start, to_print_end)
                    to_print_list.append(to_print)
                    if len(to_print_list) > 1000:
                        with open(out_file, "a") as myfile:
                            for to_print in to_print_list:
                                myfile.write(to_print)
                        to_print_list = []
                else:
                    has_next = False
            except StopIteration as e:
                has_next = False
    with open(out_file, "a") as myfile:
        for to_print in to_print_list:
            myfile.write(to_print)
    print(counter)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regions_file')
    parser.add_argument('--boundary_file')
    parser.add_argument('--out_file')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing file if exists')
    return parser.parse_args()


def main():
    args = parse_args()
    block_boundary_iterator(args.regions_file, args.boundary_file, args.out_file, args.force)


if __name__ == '__main__':
    main()
    # block_boundary_iterator("/cs/zbio/jrosensk/block_files_2/Aorta-Endothel-Z00000422.blocks.tsv",
    #                         "/cs/cbio/jon/segmentation_files/block_boundaries.tsv",
    #                         "/cs/zbio/jrosensk/block_files_2/Aorta-Endothel-Z00000422.blocks.unified.tsv")
