import datetime
import optparse
import os
from multiprocessing import Pool
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess

from asm_catalog import choose_blocks_by_fdr_bh
from asm_utils import run_command, get_tissue_file_dict


def is_bimodal(let1_uxm, let2_uxm, bimodal_thresh=0.3):
    total_u = let1_uxm[0] + let2_uxm[0]
    total_x = let1_uxm[1] + let2_uxm[1]
    total_m = let1_uxm[2] + let2_uxm[2]
    total = total_u + total_m + total_x
    prop_u = total_u / total if total > 3 else 0
    prop_m = total_m / total if total > 3 else 0

    return min(prop_u, prop_m) > bimodal_thresh


def decide_is_dirty(let1_uxm, let2_uxm, dirty_thresh=0.25):
    let1_total = let1_uxm[0] + let1_uxm[1] + let1_uxm[2]
    let2_total = let2_uxm[0] + let2_uxm[1] + let2_uxm[2]
    total_both_lets = let1_total + let2_total
    let1_prop = let1_total / total_both_lets if total_both_lets > 4 else 0
    let2_prop = let2_total / total_both_lets if total_both_lets > 4 else 0
    is_hetero = min(let1_prop, let2_prop) > 0.2
    if is_hetero:
        if is_bimodal(let1_uxm, let2_uxm):

            let1_u_prop = let1_uxm[0] / let1_total if let1_total > 3 else 0
            let1_m_prop = let1_uxm[2] / let1_total if let1_total > 3 else 0
            let2_u_prop = let2_uxm[0] / let2_total if let2_total > 3 else 0
            let2_m_prop = let2_uxm[2] / let2_total if let2_total > 3 else 0
            if min(let1_u_prop, let1_m_prop) > dirty_thresh or min(let2_u_prop, let2_m_prop) > dirty_thresh:
                return True
            else:
                return False
        else:
            return False
    else:
        return False


def decide_asm_fisher(let1_uxm, let2_uxm, psuedo_count=0, homogenous_threshold=0.75, min_read_count=2, heterozygous_prop=0):
    homogenous_relaxed = homogenous_threshold - 0.05
    counts_array = np.zeros([2, 2])
    counts_array[0, 0] = let1_uxm[0]
    counts_array[0, 1] = let1_uxm[2]
    counts_array[1, 0] = let2_uxm[0]
    counts_array[1, 1] = let2_uxm[2]
    counts_array = counts_array + psuedo_count

    oddsr, p = stats.fisher_exact(counts_array, alternative="two-sided")
    read_count_let1 = sum(let1_uxm)
    read_count_let2 = sum(let2_uxm)
    sum_um_let1 = let1_uxm[0] + let1_uxm[2]
    sum_um_let2 = let2_uxm[0] + let2_uxm[2]
    prop_let1 = read_count_let1 / (read_count_let1 + read_count_let2) if (read_count_let1 + read_count_let2) > 0 else 0
    if sum_um_let1 > min_read_count and sum_um_let2 > min_read_count and (heterozygous_prop <= prop_let1 <= (1 - heterozygous_prop)):
        u1_relaxed = let1_uxm[0] / read_count_let1
        u2_relaxed = let2_uxm[0] / read_count_let2
        m1_relaxed = let1_uxm[2] / read_count_let1
        m2_relaxed = let2_uxm[2] / read_count_let2
        u1_strict = let1_uxm[0] / sum_um_let1
        u2_strict = let2_uxm[0] / sum_um_let2
        m1_strict = let1_uxm[2] / sum_um_let1
        m2_strict = let2_uxm[2] / sum_um_let2
        homogenous1_strict = u1_strict >= homogenous_threshold and m2_strict >= homogenous_threshold and let1_uxm[
            0] > min_read_count and let2_uxm[2] > min_read_count
        homogenous1_relaxed = u1_relaxed >= homogenous_relaxed and m2_relaxed >= homogenous_relaxed
        homogenous1 = homogenous1_strict and homogenous1_relaxed
        if not homogenous1:
            homogenous2_strict = u2_strict >= homogenous_threshold and m1_strict >= homogenous_threshold and let1_uxm[
                2] > min_read_count and let2_uxm[0] > min_read_count
            homogenous2_relaxed = u2_relaxed >= homogenous_relaxed and m1_relaxed >= homogenous_relaxed
            homogenous2 = homogenous2_strict and homogenous2_relaxed
            if not homogenous2:
                p = 1
    else:
        p = -1
        homogenous1 = False
    return p, homogenous1


def decide_asm_heuristic(let1_uxm, let2_uxm, min_allele_read_count=4, homogenous_threshold=0.8,
                         similari_allele_count_thresh=0.3):
    read_count_let1 = sum(let1_uxm)
    read_count_let2 = sum(let2_uxm)
    if read_count_let1 > min_allele_read_count and read_count_let2 > min_allele_read_count:
        if (read_count_let1 > similari_allele_count_thresh * read_count_let2) and (
                read_count_let2 > similari_allele_count_thresh * read_count_let1):
            u1 = let1_uxm[0] / read_count_let1
            u2 = let2_uxm[0] / read_count_let2
            m1 = let1_uxm[2] / read_count_let1
            m2 = let2_uxm[2] / read_count_let2
            homogenous1 = u1 >= homogenous_threshold and m2 >= homogenous_threshold
            homogenous2 = m1 >= homogenous_threshold and u2 >= homogenous_threshold
            if homogenous1 or homogenous2:
                return True, homogenous1
            else:
                return False, False
    return False, False


def process_snp(snp_region, let1, let2, data_dir, homog_prop, homogenous_threshold=0.8, min_allele_read_count=4,
                similari_allele_count_thresh=0.3):
    x = 0
    chrom = snp_region.split(":")[0]
    snp_pos = int(snp_region.split(":")[1])
    bigger_region = f"{chrom}:{snp_pos - 250}-{snp_pos + 250}"
    num_times_asm_let_1_u = 0
    num_times_asm = 0
    for filename in os.listdir(data_dir):
        if filename.endswith(".bam"):
            full_file = os.path.join(data_dir, filename)
            # if not "Lung-Bronchus-Smooth-Muscle-Z00000421" in filename:
            #     continue
            # cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
            #       f"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/match_maker - |" \
            #       f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
            #       f" --min_cpg 1 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
                  f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm_3 /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
                  f" --min_cpg 1 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            res = run_command(cmd, stderr_out=subprocess.DEVNULL)

            if len(res) > 1:
                res_list = res.split("\n")
                tokens1 = res_list[0].split(": ")[1].split("\t")
                tokens2 = res_list[1].split(": ")[1].split("\t")
                tokens1 = [int(x) for x in tokens1]
                tokens2 = [int(x) for x in tokens2]
                is_asm, let1_is_u = decide_asm_heuristic(tokens1, tokens2)
                num_times_asm += int(is_asm)
                num_times_asm_let_1_u += int(let1_is_u)
    if num_times_asm >= 4 and num_times_asm_let_1_u > 0 and num_times_asm_let_1_u < num_times_asm:
        print(snp_region)
        return [snp_region]
    else:
        return []


def process_snp_fisher_with_target(snp_region, let1, let2, data_dir, homog_prop, target_tissue,
                                   homog_region_threshold=0.75):
    x = 0
    threshold_on_non_homegous = 0.1
    print(snp_region)
    chrom = snp_region.split(":")[0]
    snp_pos = int(snp_region.split(":")[1])
    bigger_region = f"{chrom}:{snp_pos - 250}-{snp_pos + 250}"
    num_times_asm_let_1_u = 0
    num_times_asm = 0
    tissue_filename_dict = get_tissue_file_dict()
    # if target_tissue in tissue_filename_dict:
    #     target_file_names = tissue_filename_dict[target_tissue]
    # else:
    #     target_file_names = []
    #     for file_list in tissue_filename_dict.values():
    #         for f in file_list:
    #             if target_tissue in f:
    #                 target_file_names.append(f)
    target_file_names = []
    for file_list in tissue_filename_dict.values():
        for f in file_list:
            if target_tissue == "Colon-Ep":
                if "Colon" in f and ("Epit" in f or "Endoc" in f):
                    target_file_names.append(f)
            elif target_tissue == "Tongue-Ep":
                if ("Ton" in f or "Thyr" in f) and "Epit" in f:
                    target_file_names.append(f)
            elif target_tissue == "Kidney-Ep":
                if "Kidney" in f and ("Epit" in f or "Pod" in f):
                    target_file_names.append(f)
            elif target_tissue == "Kidney-Ep":
                if "Small-int" in f and ("Epit" in f):
                    target_file_names.append(f)
    p_val_list = []
    target_file_names = set(target_file_names)
    bad_count = 0
    non_blank_threshold = 2
    for filename in os.listdir(data_dir):
        if filename.endswith(".bam"):
            full_file = os.path.join(data_dir, filename)
            name = filename.split(".")[0]

            # if not "Lung-Bronchus-Smooth-Muscle-Z00000421" in filename:
            #     continue
            # cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
            #       f"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/match_maker - |" \
            #       f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
            #       f" --min_cpg 1 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
                  f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
                  f" --min_cpg 3 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            res = run_command(cmd, stderr_out=subprocess.DEVNULL)

            if len(res) > 1:
                res_list = res.split("\n")
                tokens1 = res_list[0].split(": ")[1].split("\t")
                tokens2 = res_list[1].split(": ")[1].split("\t")
                tokens1 = [int(x) for x in tokens1]
                tokens2 = [int(x) for x in tokens2]
                p_val, let1_is_u = decide_asm_fisher(tokens1, tokens2, homogenous_threshold=homog_region_threshold)
                p_val_list.append((p_val, (let1_is_u, name)))
                # if name in target_file_names:
                if name in target_file_names:
                    num_blanks = get_num_blanks(non_blank_threshold, tokens1, tokens2)
                    if num_blanks != 2 and name in target_file_names:
                        bad_count += 1
    if len(p_val_list) > 0 and bad_count == 0:
        p_val_list = sorted(p_val_list, key=lambda item: item[0])
        p_val_list, let1_is_u_and_name_list = [list(t) for t in zip(*p_val_list)]

        selected_p_vals, selected_let1_is_u_and_name_list, _ = choose_blocks_by_fdr_bh(p_val_list, let1_is_u_and_name_list,
                                                                                    alpha=0.001)
        if len(selected_p_vals) > 0:
            selected_let1_is_u_list, selected_filename_list = [list(t) for t in zip(*selected_let1_is_u_and_name_list)]
            count_asm_target = 0
            for n in selected_filename_list:
                if n in target_file_names:
                    count_asm_target += 1
            num_times_asm = len(selected_let1_is_u_list)
            num_times_asm_let_1_u = sum(selected_let1_is_u_list)
            if num_times_asm >= 2 and 0 < num_times_asm_let_1_u < num_times_asm:
                print(f"imp\t{snp_region}")
                return [snp_region]
            else:
                return []
        else:
            return []
    else:
        return []


def get_num_blanks(non_blank_threshold, tokens1, tokens2):
    num_blanks = 4
    if tokens1[0] > non_blank_threshold:
        num_blanks -= 1
    if tokens1[2] > non_blank_threshold:
        num_blanks -= 1
    if tokens2[0] > non_blank_threshold:
        num_blanks -= 1
    if tokens2[2] > non_blank_threshold:
        num_blanks -= 1
    return num_blanks


def classify_as_imprinting(num_times_asm, num_times_asm_let_1_u):
    return num_times_asm >= 3 and 0 < num_times_asm_let_1_u < num_times_asm


def classify_as_sd(num_times_asm, num_times_asm_let_1_u, count_bad):
    return num_times_asm >= 3 and (0 == num_times_asm_let_1_u or num_times_asm_let_1_u == num_times_asm) and (
                count_bad < 20)


def process_snp_fisher(snp_region, let1, let2, data_dir, homog_prop, homog_region_threshold=0.75, find_imp=True):
    x = 0
    print(snp_region)
    chrom = snp_region.split(":")[0]
    snp_pos = int(snp_region.split(":")[1])
    bigger_region = f"{chrom}:{snp_pos - 250}-{snp_pos + 250}"
    num_times_asm_let_1_u = 0
    num_times_asm = 0
    non_blank_threshold = 2
    p_val_list = []
    count_bad_sd = 0
    count_is_dirty = 0
    count_bimodal = 0
    asm_files = []
    for filename in os.listdir(data_dir):
        if filename.endswith(".bam"):
            full_file = os.path.join(data_dir, filename)
            # if not "Lung-Bronchus-Smooth-Muscle-Z00000421" in filename:
            #     continue
            # cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
            #       f"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/match_maker - |" \
            #       f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
            #       f" --min_cpg 1 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
                  f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/CpG.blacklist_removed.bed.gz {bigger_region}" \
                  f" --min_cpg 3 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            res = run_command(cmd, stderr_out=subprocess.DEVNULL)

            if len(res) > 1:
                res_list = res.split("\n")
                tokens1 = res_list[0].split(": ")[1].split("\t")
                tokens2 = res_list[1].split(": ")[1].split("\t")
                tokens1 = [int(x) for x in tokens1]
                tokens2 = [int(x) for x in tokens2]
                p_val, let1_is_u = decide_asm_fisher(tokens1, tokens2, homogenous_threshold=homog_region_threshold)
                num_blanks = get_num_blanks(non_blank_threshold, tokens1, tokens2)
                is_good = num_blanks == 2 and (
                            (tokens1[0] >= 3 and tokens2[2] >= 3) or (tokens1[2] >= 3 and tokens2[0] >= 3))

                # is_good_sd = num_blanks >= 2 and (
                #             (tokens1[0] >= 5 and tokens1[2] <= 3) or (tokens1[0] <= 3 and tokens1[2] >= 5) or (
                #                 tokens2[0] >= 5 and tokens2[2] <= 3) or (tokens2[0] <= 3 and tokens2[2] >= 5) or (
                #                         sum(tokens1) < 3 and sum(tokens2) < 3))
                # count_bad_sd += int(not is_good_sd)
                is_bimodal_a = is_bimodal(tokens1, tokens2)
                count_bimodal += is_bimodal_a
                is_dirty = decide_is_dirty(tokens1, tokens2)
                count_is_dirty += int(is_dirty)
                p_val_list.append((p_val, (let1_is_u, (is_good, filename.split(".")[0]))))
    if len(p_val_list) > 0:
        p_val_list = sorted(p_val_list, key=lambda item: item[0])
        p_val_list = [(p_val, (a, (is_2, f_name))) for p_val, (a, (is_2, f_name)) in p_val_list if (is_2 and p_val < 1)]
        if len(p_val_list) > 0:
            p_val_list, let1_is_u_list_and_2_blanks = [list(t) for t in zip(*p_val_list)]
            let1_is_u_list, is_2_blanks_listand_filenames = [list(t) for t in zip(*let1_is_u_list_and_2_blanks)]
            is_2_blanks, filenames_sig = [list(t) for t in zip(*is_2_blanks_listand_filenames)]
            let1_u_list_and_fnames = list(zip(let1_is_u_list, filenames_sig))
            selected_p_vals, selected_let1_is_u_list_and_f_names = p_val_list, let1_u_list_and_fnames
            # selected_p_vals, selected_let1_is_u_list_and_f_names = choose_blocks_by_fdr_bh(p_val_list,
            #                                                                                let1_u_list_and_fnames,
            #                                                                                alpha=0.1)
            if len(selected_p_vals) == 0:
                return []
            selected_let1_is_u_list, filenames_sig = [list(t) for t in zip(*selected_let1_is_u_list_and_f_names)]
            num_times_asm = len(selected_let1_is_u_list)
            num_times_asm_let_1_u = sum(selected_let1_is_u_list)
            if find_imp:
                to_print = f"{snp_region}\t{','.join(filenames_sig)}\t{','.join([str(p) for p in selected_p_vals])}\t" \
                           f"{','.join([str(l_is_u) for l_is_u in selected_let1_is_u_list])}\t{count_is_dirty}\t{count_bimodal}"
                if classify_as_imprinting(num_times_asm, num_times_asm_let_1_u):

                    print(f"imp\t{to_print}")
                    return [to_print]
                else:
                    return []
                # elif classify_as_sd(num_times_asm, num_times_asm_let_1_u, count_bad_sd):
                #     print(f"sd-asm\t{to_print}")
                #     return []
                # else:
                #     print(f"snp\t{to_print}")
                #     return []
            else:
                if classify_as_sd(num_times_asm, num_times_asm_let_1_u, count_bad_sd):
                    print(f"sd-asm\t{snp_region}")
                    return [snp_region]
                else:
                    return []
            #
            # if num_times_asm >= 2 and 0 < num_times_asm_let_1_u < num_times_asm:
            #     print(f"imp\t{snp_region}")
            #     return [snp_region]
            # else:
            #     return []
        else:
            return []
    else:
        return []


def print_snp(snp_region, let1, let2, data_dir, homog_prop, num_treahds=16):
    chrom = snp_region.split(":")[0]
    snp_pos = int(snp_region.split(":")[1])
    bigger_region = f"{chrom}:{snp_pos - 250}-{snp_pos + 250}"
    num_times_asm_let_1_u = 0
    num_times_asm = 0
    p_val_list = []
    params = []
    for filename in os.listdir(data_dir):
        if filename.endswith(".bam"):
            full_file = os.path.join(data_dir, filename)
            # if not "Ovary-Epithelial-Z000000QT" in filename:
            #     continue
            # cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
            #       f"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/match_maker - |" \
            #       f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
            #       f" --min_cpg 1 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            cmd = f"samtools view {full_file} {bigger_region} -F 1796 -q 10 -f 3 | " \
                  f" /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/pipeline_wgbs/summary_uxm /cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/CpG.bed.gz {bigger_region}" \
                  f" --min_cpg 3 --clip 0 --snp_pos {snp_pos} --snp_let1 '{let1}' --snp_let2 '{let2}' --homog_prop {homog_prop}"
            params.append((cmd, filename, snp_region, let1, let2, homog_prop))

    with Pool(num_treahds) as p:
        res = p.starmap(handle_single_cmd, params)
        p.close()
        p.join()

    for to_print in res:
        print(to_print)


def handle_single_cmd(cmd, filename, snp_region, snp_let1, snp_let2, homog_thresh=0.65):
    res = run_command(cmd, stderr_out=subprocess.DEVNULL)
    if len(res) > 1:
        res_list = res.split("\n")
        tokens1 = res_list[0].split(": ")[1].split("\t")
        tokens2 = res_list[1].split(": ")[1].split("\t")
        tokens1 = [int(x) for x in tokens1]
        tokens2 = [int(x) for x in tokens2]
        total1 = sum(tokens1)
        uprop1 = tokens1[0] / total1 if total1 > 0 else 0
        mprop1 = tokens1[2] / total1 if total1 > 0 else 0
        total2 = sum(tokens2)
        uprop2 = tokens2[0] / total2 if total2 > 0 else 0
        mprop2 = tokens2[2] / total2 if total2 > 0 else 0
        p_val, let1_is_u = decide_asm_fisher(tokens1, tokens2, homogenous_threshold=0.65)
        amount_str1 = f"{tokens1[0]}/{tokens1[2]}"
        if total1 > 4:
            if mprop1 >= homog_thresh:
                amount_str1 = f"\033[1;31m{tokens1[0]}/{tokens1[2]}\033[0m"
            elif uprop1 >= homog_thresh:
                amount_str1 = f"\033[1;32m{tokens1[0]}/{tokens1[2]}\033[0m"
            else:
                amount_str1 = f"\033[1;93m{tokens1[0]}/{tokens1[2]}\033[0m"
        amount_str2 = f"{tokens2[0]}/{tokens2[2]}"
        if total2 > 4:
            if mprop2 >= homog_thresh:
                amount_str2 = f"\033[1;31m{tokens2[0]}/{tokens2[2]}\033[0m"
            elif uprop2 >= homog_thresh:
                amount_str2 = f"\033[1;32m{tokens2[0]}/{tokens2[2]}\033[0m"
            else:
                amount_str2 = f"\033[1;93m{tokens2[0]}/{tokens2[2]}\033[0m"
        res_str = f"{filename}\t{snp_region}\t{amount_str1}\t{amount_str2}\t{p_val}\t{snp_let1}\t{snp_let2}"
        return res_str
    return ""


if __name__ == '__main__':
    # counts_array = np.zeros([2, 2])
    # counts_array[0, 0] = 46
    # counts_array[0, 1] = 0
    # counts_array[1, 0] = 0
    # counts_array[1, 1] = 30
    # # counts_array = counts_array + psuedo_count
    # 
    # oddsr, p = stats.fisher_exact(counts_array, alternative="two-sided")
    parser = optparse.OptionParser()
    parser.add_option('--snp_pos',
                      default="")
    parser.add_option('--snp_let_1',
                      default="T")
    parser.add_option('--snp_let_2',
                      default="G")
    parser.add_option('--snp_file',
                      default="")
    parser.add_option('--num_threads',
                      default=32)
    parser.add_option('--homog_prop',
                      default=0.65)
    parser.add_option('--target_tissue',
                      default="")
    parser.add_option('--find_sd', action="store_true")
    options, arguments = parser.parse_args()
    num_treahds = int(options.num_threads)
    # tabix -R /cs/zbio/jrosensk/block_files_smoothed/only_colon/bimodal_in_all_colon.not_in_papers.bed gnomAD.all.bed.gz | awk '($7 > 0.1 && $7 < 0.9)' > /cs/zbio/jrosensk/block_files_smoothed/only_colon/snp_in_bimodal_colon_blocks.bed
    print("{}: starting handle snps".format(str(datetime.datetime.now())))
    # snp_file = "/cs/zbio/jrosensk/block_files_smoothed/only_tongue_epithelial/snp_in_bimodal_tongue_ep_blocks.bed"
    # snp_file = "/cs/zbio/jrosensk/tmp_imprinting_files/snps_near_imp.not_near_known_icr.2.bed"
    # snp_file = "/cs/zbio/jrosensk/tmp_imprinting_files/small_2.gnom.bed"
    snp_file = options.snp_file
    target_tissue = options.target_tissue
    params = []
    # with open(snp_file, 'r') as f:
    #     for line in f:
    #         line = line.strip()
    #         tokens = line.split("\t")
    #         snp_pos = f"{tokens[0]}:{tokens[1]}"
    #         params.append((snp_pos, tokens[3], tokens[4], "/cs/cbio/jon/grail_atlas/data", 0.65))
    # process_snp_fisher("chr7:76253066", 'C', 'T', "/cs/cbio/jon/grail_atlas/data", 0.65)
    # process_snp_fisher("chr11:2720873", 'C', 'G', "/cs/cbio/jon/grail_atlas/data", 0.65)
    find_imp = not options.find_sd
    if options.snp_file != "" and options.snp_pos == "":
        with open(snp_file, 'r') as f:
            for line in f:
                line = line.strip()
                tokens = line.split("\t")
                snp_pos = f"{tokens[0]}:{tokens[1]}"
                # snp_pos = "chr15:23811592"
                # tokens[3] = 'C'
                # tokens[4] = 'T'
                params.append((snp_pos, tokens[3], tokens[4], "/cs/cbio/jon/grail_atlas/data", 0.65, 0.75, find_imp))

        with Pool(num_treahds) as p:
            res = p.starmap(process_snp_fisher, params)
            p.close()
            p.join()
        # with Pool(num_treahds) as p:
        #     res = p.starmap(process_snp_fisher_with_target, params)
        #     p.close()
        #     p.join()
        # # # # with Pool(num_treahds) as p:
        # # # #     res = p.starmap(process_snp_fisher, params)
        # # # #     p.close()
        # # # #     p.join()
        imp_snp_list = [item for sublist in res for item in sublist]
        print("printing all imp snps")
        for snp in imp_snp_list:
            print(snp)

    # process_snp("chr12:31271788", 'G', 'A', "/cs/cbio/jon/grail_atlas/data", 0.65)

    if options.snp_file == "" and options.snp_pos != "":
        print_snp(options.snp_pos, options.snp_let_1, options.snp_let_2, "/cs/cbio/jon/grail_atlas/data", float(options.homog_prop),
                  num_treahds=num_treahds)
    # print_snp("chr11:2017110", 'C', 'G', "/cs/cbio/jon/grail_atlas/data", 0.65,
    #           num_treahds=1)
    # process_snp_fisher("chr21:15352485", 'G', 'T', "/cs/cbio/jon/grail_atlas/data", 0.65)
    # process_snp_fisher("chr2:81162067", 'C', 'T', "/cs/cbio/jon/grail_atlas/data", 0.65)
    # process_snp_fisher("chr2:81162072", 'A', 'C', "/cs/cbio/jon/grail_atlas/data", 0.65)
    # process_snp_fisher("chr2:81409135", 'T', 'A', "/cs/cbio/jon/grail_atlas/data", 0.65)
    # process_snp_fisher("chr1:17019022", 'T', 'C', "/cs/cbio/jon/grail_atlas/data", 0.65)
    x = 0
    print("{}: finished handle snps".format(str(datetime.datetime.now())))
