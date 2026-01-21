import pandas as pd


if __name__ == '__main__':
    
    blat_regions_file = "sd_asm_analysis/homog/homog_aligned/all_snps/meqtl/simple_sd_asm/blat_found.txt"

    regions_dict = {}

    with open(blat_regions_file, 'r') as f:
        for line in f:
            reg = line.strip()
            if reg in regions_dict:
                cur_count = regions_dict[reg]
            else:
                cur_count = 0
            cur_count += 1
            regions_dict[reg] = cur_count

    multi_cnt_regions = [re for re, cnt in list(regions_dict.items()) if cnt > 1]
    not_multi_regions = [re for re, cnt in list(regions_dict.items()) if cnt == 1]

    with open("sd_asm_analysis/homog/homog_aligned/all_snps/meqtl/simple_sd_asm/blat_duplicate_regions.txt", 'w') as f_out:
        for rg in multi_cnt_regions:
            chrom = rg.split(":")[0]
            s_e = rg.split(":")[1]
            start = s_e.split("-")[0]
            end = s_e.split("-")[1]
            f_out.write(f"{chrom}\t{start}\t{end}\n")

    x = 0
