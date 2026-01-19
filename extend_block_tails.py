"""
Bimodal Methylation Block Boundary Refinement Tool
------------------------------------------
This script refines the start and end boundaries (CpG indices) of previously
identified bimodal DNA methylation blocks. It "explores the tails" of each block by
analyzing raw fragment patterns (PAT files) and methylation levels (BETA files)
to extend or trim boundaries based on bimodal signals.

Key Capabilities:
- Tail Exploration: Heuristically adjusts block borders by scanning adjacent
  CpG sites for persistent bimodal (U/M) or heterogeneous (X) read proportions.
- Tabix Integration: Uses indexed genomic queries to efficiently fetch local
  methylation fragments.
- Multiprocessing: Processes genomic regions in parallel across chromosomes
  to handle large-scale WGBS data.
- Signal Validation: Implements error-tolerant scanning logic (allowed mistakes)
  to prevent premature boundary termination due to single-site noise.
"""

import optparse
import os
import subprocess
import warnings
from multiprocessing import Pool

import numpy as np
import math

from asm_utils import get_chroms


def explore_tail_bimodal_finder(chrom, start_site, end_site, pat_file, beta_file):
    max_addition = 25
    accepted_threshold = 0.35
    allowed_mistakes = 5
    new_block_min_len = 3
    allowed_x_prop = 0.4
    min_um_proportion = 0.3
    homog_prop = 0.25
    new_found_end = find_border_by_um_prop(allowed_mistakes, min_um_proportion, chrom, end_site, end_site + max_addition, homog_prop,
                                           pat_file, 1)
    new_found_start = find_border_by_x_prop(allowed_mistakes, min_um_proportion, chrom, start_site,
                                            start_site - max_addition, homog_prop,
                                            pat_file, -1)

    return min(new_found_start, start_site), max(new_found_end, end_site) # this is necessary because sometimes the very first site breaks rules and is removed


def explore_tail(chrom, start_site, end_site, pat_file, beta_file):
    max_addition = 25
    accepted_threshold = 0.35
    allowed_mistakes = 1
    new_block_min_len = 3
    allowed_x_prop = 0.4
    homog_prop = 0.25

    beta_array = np.fromfile(beta_file, dtype=np.uint8).reshape((-1, 2))
    found_end = find_new_border_by_c_prop(accepted_threshold, allowed_mistakes, beta_array, end_site, max_addition, 1)
    if found_end - end_site >= new_block_min_len:
        new_found_end = find_border_by_tabix(chrom, end_site, found_end, pat_file, 1)
        # new_found_end = found_end
    else:
        new_found_end = end_site
    if new_found_end == -1:
        print(f"{start_site}-{end_site}")
        new_found_end = end_site

    found_start = find_new_border_by_c_prop(accepted_threshold, allowed_mistakes, beta_array, start_site, max_addition, -1)
    if start_site - found_start >= new_block_min_len:
        # new_found_start = find_border_by_x_prop(allowed_mistakes, allowed_x_prop, chrom, start_site, found_start, homog_prop,
        #                                       pat_file, -1)
        new_found_start = find_border_by_tabix(chrom, start_site, found_start, pat_file, -1)
        # new_found_start = found_start
    else:
        new_found_start = start_site
    if new_found_start == -1:
        print(f"{start_site}-{end_site}")
        new_found_start = end_site

    return min(new_found_start, start_site), max(new_found_end, end_site) # this is necessary because sometimes the very first site breaks rules and is removed


def find_border_by_um_prop(allowed_mistakes, allowed_min_um_prop, chrom, border_site, border_extension, homog_prop, pat_file, direction):
    extension_length = (direction * border_extension) - (border_site * direction)
    u_count_array = np.zeros((extension_length))
    m_count_array = np.zeros((extension_length))
    total_count_array = np.zeros((extension_length))
    region_to_pull = f"{chrom}:{border_site}-{border_extension}" if direction == 1 else f"{chrom}:{border_extension}-{border_site}"
    tabix_cmd = f"tabix {pat_file} {region_to_pull}"
    lines = subprocess.check_output(tabix_cmd, shell=True).decode().split("\n")
    beginning = border_site if direction == 1 else border_extension
    is_first = True
    for line in lines:
        if len(line) < 2:
            break
        tokens = line.split("\t")
        if is_first:
            first_site = int(tokens[1])
            is_first = False
        last_site = int(tokens[1])
        pat = tokens[2]
        c_count = pat.count("C")
        t_count = pat.count("T")
        total = (c_count + t_count)
        c_prop = c_count / total if total > 0 else 0
        if total > 0:
            cur_start = last_site - beginning
            cur_end = min(extension_length, (cur_start + len(pat)))
            total_count_array[cur_start:cur_end] += 1
            # if cur_start in range(5, 10):
            #     x = 0
            if c_prop <= homog_prop:
                u_count_array[cur_start:cur_end] += 1
            elif c_prop >= (1 - homog_prop):
                m_count_array[cur_start:cur_end] += 1
            # if homog_prop < c_prop < (1 - homog_prop):
            #     x_count_array[cur_start:cur_end] += 1
    num_mistakes = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u_props = u_count_array / total_count_array
        m_props = m_count_array / total_count_array
    u_props[np.isinf(u_props)] = 0
    u_props[np.isnan(u_props)] = 0
    m_props[np.isinf(m_props)] = 0
    m_props[np.isnan(m_props)] = 0
    if not is_first:
        border_extension = min(border_extension, last_site) if direction == 1 else max(border_extension, first_site)
        new_extension_length = (direction * border_extension) - (border_site * direction) # to handle cross chromosome sites
    else:
        new_extension_length = extension_length
    if new_extension_length != extension_length:
        diff = extension_length - new_extension_length
        if direction == 1:
            u_props = u_props[:-diff]
            m_props = m_props[:-diff]
        else:
            u_props = u_props[diff:]
            m_props = m_props[diff:]
    directional_range = range(0, new_extension_length) if direction == 1 else range(new_extension_length-1, -1, direction)
    num_steps = 0
    for i in directional_range:
        cur_prop = min(u_props[i], m_props[i])
        if cur_prop < allowed_min_um_prop:
            num_mistakes += 1
        else:
            num_mistakes = 0
        if num_mistakes > allowed_mistakes:
            break
        num_steps += 1
    new_found_end = border_site + (direction * (num_steps + 1)) - (num_mistakes * direction)
    return new_found_end


def find_border_by_tabix(chrom, border_site, border_extension, pat_file, direction, max_read_len=50):
    region_to_pull = f"{chrom}:{border_site}-{border_extension}" if direction == 1 else f"{chrom}:{border_extension - max_read_len}-{border_site}"
    tabix_cmd = f"tabix {pat_file} {region_to_pull}"
    lines = subprocess.check_output(tabix_cmd, shell=True).decode().split("\n")
    is_first = True
    last_site = -1
    first_site = -1
    for line in lines:
        if len(line) < 2:
            break
        tokens = line.split("\t")
        if is_first:
            first_site = int(tokens[1])
            is_first = False
        last_site = int(tokens[1])

    if direction == 1:
        new_found_end = min(border_extension, last_site)
    else:
        new_found_end = max(border_extension, first_site)
    return new_found_end


def find_border_by_x_prop(allowed_mistakes, allowed_x_prop, chrom, border_site, border_extension, homog_prop, pat_file, direction, max_read_len=50):
    extension_length = (direction * border_extension) - (border_site * direction)
    x_count_array = np.zeros((extension_length))
    total_count_array = np.zeros((extension_length))
    region_to_pull = f"{chrom}:{border_site}-{border_extension}" if direction == 1 else f"{chrom}:{border_extension-max_read_len}-{border_site}"
    tabix_cmd = f"tabix {pat_file} {region_to_pull}"
    lines = subprocess.check_output(tabix_cmd, shell=True).decode().split("\n")
    beginning = border_site if direction == 1 else border_extension
    is_first = True
    last_site = -1
    first_site = -1
    for line in lines:
        if len(line) < 2:
            break
        tokens = line.split("\t")
        if is_first:
            first_site = int(tokens[1])
            is_first = False
        last_site = int(tokens[1])
        pat = tokens[2]
        c_count = pat.count("C")
        t_count = pat.count("T")
        total = (c_count + t_count)
        c_prop = c_count / total if total > 0 else 0
        if total > 0:
            cur_start = last_site - beginning
            cur_end = min(extension_length, (cur_start + len(pat)))
            if direction == -1 and cur_end < border_extension:
                continue
            total_count_array[cur_start:cur_end] += 1
            if homog_prop < c_prop < (1 - homog_prop):
                x_count_array[cur_start:cur_end] += 1
    num_mistakes = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x_props = x_count_array / total_count_array
    x_props[np.isinf(x_props)] = 0
    x_props[np.isnan(x_props)] = 0
    if not is_first:
        border_extension = min(border_extension, last_site) if direction == 1 else max(border_extension, first_site)
        new_extension_length = (direction * border_extension) - (border_site * direction) # to handle cross chromosome sites
    else:
        new_extension_length = extension_length
    if new_extension_length != extension_length:
        diff = extension_length - new_extension_length
        if direction == 1:
            x_props = x_props[:-diff]
        else:
            x_props = x_props[diff:]
    directional_range = range(0, new_extension_length) if direction == 1 else range(new_extension_length-1, -1, direction)
    num_steps = 0
    for i in directional_range:
        cur_prop = x_props[i]
        if cur_prop > allowed_x_prop:
            num_mistakes += 1
        else:
            num_mistakes = 0
        if num_mistakes > allowed_mistakes:
            break
        num_steps += 1
    new_found_end = border_site + (direction * (num_steps+1)) - (num_mistakes * direction)
    if last_site == -1 or first_site == -1:
        return -1
    if direction == 1:
        new_found_end = min(new_found_end, last_site)
    else:
        new_found_end = max(new_found_end, first_site)
    return new_found_end


def find_new_border_by_c_prop(accepted_threshold, allowed_mistakes, beta_array, border_site, max_addition, direction):
    found_end = border_site
    num_mistakes = 0
    consecutive_mistakes = 0
    consecutive_success = 0
    range_end = found_end + (max_addition * direction)
    for cur_ind in range(found_end, range_end, direction):
        row = beta_array[cur_ind]
        methyl_prop = row[0] / row[1] if row[1] > 0 else 0
        if not math.isnan(methyl_prop):
            if not accepted_threshold < methyl_prop < (1 - accepted_threshold):
                num_mistakes += 1
                consecutive_mistakes += 1
                consecutive_success = 0
            else:
                consecutive_success += 1
                consecutive_mistakes = 0
            if consecutive_success > 10:
                num_mistakes = 0
        if consecutive_mistakes > allowed_mistakes:
            break
        found_end = found_end + (1 * direction)
    found_end = found_end - (consecutive_mistakes * direction)
    return found_end


def process_region_blocks_file(region, blocks_file, pat_file, beta_file):
    tabix_cmd = f"tabix {blocks_file} {region}"
    cur_blocks_lines = subprocess.check_output(tabix_cmd, shell=True).decode().split("\n")
    new_regions = []
    for line in cur_blocks_lines:
        if len(line) > 2:
            tokens = line.split("\t")
            chrom = tokens[0]
            start_site = int(tokens[3])
            end_site = int(tokens[4])
            new_start_site, new_end_site = explore_tail(chrom, start_site, end_site, pat_file, beta_file)
            new_regions.append((new_start_site, new_end_site))
    print(f"finished processing {region}")
    return new_regions


def process_blocks_file(blocks_file, pat_file, beta_file, out_file, num_threads):
    chroms_list = get_chroms()
    params_list = [(chrom, blocks_file, pat_file, beta_file) for chrom in chroms_list]
    with Pool(num_threads) as p:
        res = p.starmap(process_region_blocks_file, params_list)
        p.close()
        p.join()
    new_regions = [item for sublist in res for item in sublist]
    with open(out_file, "w") as f_out:
        for regions_sites in new_regions:
            f_out.write(f"{regions_sites[0]}\t{regions_sites[1]}\n")


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--in_file', require=True)
    parser.add_option('--data_dir', require=True)
    parser.add_option('--num_threads',
                      default=1)
    parser.add_option('--out_file', require=True)
    options, arguments = parser.parse_args()
    file_name = options.in_file

    base_name = file_name.split("/")[-1].split(".")[0]
    pat_file = os.path.join(options.data_dir, base_name + ".pat.gz")
    beta_file = os.path.join(options.data_dir, base_name + ".beta")
    process_blocks_file(file_name, pat_file, beta_file, options.out_file, int(options.num_threads))