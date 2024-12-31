import optparse

import numpy as np
import pandas as pd
import scipy.stats as stats
import gzip as gz
from handle_snps import decide_asm_fisher, get_num_blanks
import os

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('--in_file',
                      default="")
    parser.add_option('--data_dir',
                      default="")
    parser.add_option('--file_suffix',
                      default=".allele.homog.gz")
    options, arguments = parser.parse_args()
    file_name = options.in_file
    base_name = os.path.basename(file_name).split(".")[0]
    sample_name = base_name.split(".")[0]
    # input_dir = "/cs/zbbio/jrosensk/blocks_masked_imprinting/homog/homog_aligned/all_snps"
    file_suffix = options.file_suffix
    expected_cols = ["chrom", "start", "let1", "let2", "u1", "x1", "m1", "u2", "x2", "m2"]
    # out_path = os.path.join(input_dir, f"{sample_name}{file_suffix}")
    # for filename in os.listdir(input_dir):
    #     if not filename.endswith(file_suffix):
    #         continue
    #     full_path = os.path.join(input_dir, filename)
    non_blank_threshold = 2
    header_df = pd.read_csv(file_name, sep="\t", header=None, nrows=3)
    if len(header_df.columns) > len(expected_cols):
        print(f"{sample_name} has more columns than expected")
        raise Exception
    p_val_list = []
    is_let1_u_list = []
    with gz.open(file_name, 'r') as f:
        for line in f:
            line = line.decode().strip()
            tokens = line.split("\t")
            snp = f"{tokens[0]}:{tokens[1]}"
            # if snp != "chr1:39559408":
            #     continue
            uxm1 = [int(x) for x in tokens[4:7]]
            uxm2 = [int(x) for x in tokens[7:]]
            p_val, let1_is_u = decide_asm_fisher(uxm1, uxm2, homogenous_threshold=0.75, heterozygous_prop=0.1)
            num_blanks = get_num_blanks(non_blank_threshold, uxm1, uxm2)
            is_good = num_blanks == 2 and (
                    (uxm1[0] >= 3 and uxm2[2] >= 3) or (uxm1[2] >= 3 and uxm2[0] >= 3))
            if is_good:
                p_val_list.append(p_val)
                is_let1_u_list.append(let1_is_u)
            else:
                p_val_list.append(-1)
                is_let1_u_list.append(False)

    cur_df = pd.read_csv(file_name, sep="\t", header=None,
                         names=expected_cols)
    cur_df['p_val'] = p_val_list
    cur_df['let1_is_u'] = is_let1_u_list
    out_path = os.path.join(options.data_dir, f"{base_name}.allele.homog.pval.gz")
    cur_df.to_csv(out_path, sep="\t", index=False, header=False)
