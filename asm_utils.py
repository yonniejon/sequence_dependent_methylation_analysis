import pickle
import os.path
from enum import Enum
import subprocess
import pandas as pd
import shutil
from pathlib import Path
import uuid

imprinted_genes_with_placenta = set(
        ["DIRAS3", "IL12RB2", "RNU5D-1", "AGO1", "UTS2", "THAP3", "CACNA1E", "CYP2J2", "ACOT11", "LINC00467", "LRRTM1",
         "TMEM247", "THUMPD2", "PAX8", "PAX8-AS1", "DNAH7", "ICA1L", "GPR1", "GPR1-AS", "ZDBF2", "MRPL44", "SPHKAP",
         "USP4", "SLC4A7", "ZNF385D", "EFCC1", "RAB7A", "MCCC1", "FGF12", "MFI2", "GPR78", "STX18-AS1", "PDE6B",
         "SH3BP2", "NAP1L5", "GRID2", "SFRP2", "FAM149A", "FRG1", "PLEKGH4B", "RHOBTB3", "NUDT12", "VTRNA2-1",
         "ZNF354C", "CUL7", "MDGA1", "MOCS1", "C6orf47", "RNF144B", "CD83", "FAM50B", "FAM50B-AS", "AIM1", "LIN28B",
         "PHACTR2", "HYMAI", "PLAGL1", "SLC22A2", "SLC22A3", "PLG", "KIF25", "GRB10", "RAPGEF5", "SCIN", "THSD7A",
         "CALCR", "TFPI2", "SGCE", "PEG10", "PDK4", "CPA4", "MEST", "MESTIT1", "COPG2IT1", "KLF14", "KLHDC10", "AGBL3",
         "PRKAG2", "PTK2B", "R3HCC1", "CLDN23", "DLGAP2", "PKIA", "ZFAT", "ZFAT-AS1", "PEG13", "PSCA", "NAPRT1",
         "TRAPPC9", "KCNK9", "DENND3", "GLIS3", "PGM5P3-AS1", "EXD3", "PTCHD3", "ITGA8", "PROSER2", "PROSER2-AS1",
         "JMJD1C", "AIFM2", "USMG5", "VWA2", "INPP5F_V2", "CPXM2", "ACCS", "ALKBH3", "MAPK8IP1", "WT1-Alttranscript",
         "WT1AS", "LINC00294", "miR-675", "IGF2", "IGF2AS", "INS", "KCNQ1", "KCNQ1OT1", "KCNQ1DN", "CDKN1C", "PHLDA2",
         "miR-483", "SLC22A18", "ZNF215", "NAV2", "ART5", "OVCH2", "RNF141", "IRF7", "ANO1", "PAK1", "VSTM5", "ZC3H12C",
         "SPA17", "NTM", "OPCML", "TIGAR", "CACNA1C", "WIF1", "N4BP2L1", "RB1", "RB2", "LPAR6", "DLEU7", "KLHL1",
         "FGF14", "PCK2", "PAPLN-AS1", "DLK1", "MEG3", "MIR337", "RTL1", "MEG8", "miR-134", "PiRNAs", "MKRN3", "MAGEL2",
         "NDN", "NPAP1", "SNURF", "SNRPN", "SNORD107", "SNORD64", "SNORD108", "SNORD109A", "SNORD116@", "IPW",
         "SNORD115@", "SNORD109B", "UBE3A-AS", "UBE3A", "PWRN1", "SNHG14", "H73492", "RYR3", "DNM1P35", "RASGRF1",
         "FAM174B", "IRAIN", "LRRK1", "SIAH1", "ZNF597", "NAA60isoform1", "PDPR", "ZFP90", "CLEC3A", "NLGN2", "SEPT4",
         "ZNF714", "AXL", "DNMT1", "SIPR2", "ICAM1", "FDX1L", "ZNF833P", "GNG7", "ANO8", "CACNA1A", "C19MC", "ZNF331",
         "Anti-MIR371-MIR373", "MIR512-1", "ZIM2", "PEG3", "MIMT1", "ZNF542P", "CST1", "PSIMCT-1", "ACTL10", "NNAT",
         "BLCAP", "ZHX3", "L3MBTL", "SGK2", "CYP24A1", "NESP55", "GNASXL", "Exon-1A", "GS-alpha", "SANG", "miR-296",
         "miR-298", "GDAP1L1", "PRMT2", "CBR1", "TPTEP1", "ARVCF", "CACNA1l", "NHP2L1", "SLC9A7"])

imprinted_genes = set(['DIRAS3', 'RNU5D-1', 'UTS2', 'LRRTM1', 'PAX8', 'PAX8-AS1', 'GPR1', 'ZDBF2', 'FRG1', 'RHOBTB3', 'VTRNA2-1', 'FAM50B', 'FAM50B-AS', 'PHACTR2', 'HYMAI', 'PLAGL1', 'SLC22A2', 'SLC22A3', 'KIF25', 'GRB10', 'CALCR', 'TFPI2', 'SGCE', 'PEG10', 'CPA4', 'MEST', 'MESTIT1', 'COPG2IT1', 'KLF14', 'ZFAT-AS1', 'PEG13', 'NAPRT1', 'TRAPPC9', 'KCNK9', 'INPP5F_V2', 'WT1-Alttranscript', 'WT1AS', 'IGF2', 'IGF2AS', 'INS', 'KCNQ1', 'KCNQ1OT1', 'KCNQ1DN', 'CDKN1C', 'PHLDA2', 'miR-483', 'SLC22A18', 'ZNF215', 'ANO1', 'NTM', 'WIF1', 'RB1', 'RB2', 'LPAR6', 'DLK1', 'MEG3', 'MIR337', 'RTL1', 'MEG8', 'miR-134', 'PiRNAs', 'MKRN3', 'MAGEL2', 'NDN', 'NPAP1', 'SNURF', 'SNRPN', 'SNORD107', 'SNORD64', 'SNORD108', 'SNORD109A', 'SNORD116@', 'IPW', 'SNORD115@', 'SNORD109B', 'UBE3A-AS', 'UBE3A', 'PWRN1', 'SNHG14', 'H73492', 'IRAIN', 'ZNF597', 'NAA60isoform1', 'AXL', 'GNG7', 'C19MC', 'ZNF331', 'Anti-MIR371-MIR373', 'ZIM2', 'PEG3', 'MIMT1', 'ZNF542P', 'CST1', 'PSIMCT-1', 'NNAT', 'BLCAP', 'L3MBTL', 'SGK2', 'NESP55', 'GNASXL', 'Exon-1A', 'GS-alpha', 'SANG', 'miR-296', 'miR-298', 'GDAP1L1', 'SLC9A7', 'GNAS','NAA60','DLK-1','DIO-3', "H19", "AIRN", "DLK-1", "GNAS"])
imprinted_genes = set([gene.upper() for gene in imprinted_genes])


def delete_dir(dir_name):
    shutil.rmtree(dir_name)

def create_unique_dir_name():
    cur_id = uuid.uuid1()
    path_name = str(cur_id) + "_tmp"

    Path(path_name).mkdir(parents=True, exist_ok=False)
    return path_name


def create_dir(name):
    path_name = name + "_tmp"

    Path(path_name).mkdir(parents=True, exist_ok=True)
    return path_name


def get_sample_name(file_name, suffix=".uxm.bed.gz"):
    basename = os.path.basename(file_name)
    if file_name.endswith(suffix):
        sample_name = basename[:-len(suffix)]
    return sample_name


def open_pickle_object(file_name):
    with open(file_name, 'rb') as handle:
        res = pickle.load(handle)
    return res


def format_float(value, num_digits=3):
    if 0 < value < 1 and value < 0.001:
        format_string = "{:." + str(num_digits) + "e}"
        return format_string.format(value)
    else:
        format_string = "{:." + str(num_digits) + "f}"
        return format_string.format(value)


def format_float_percent(value, num_digits=3):
    if 0 < value < 1 and value < 0.001:
        format_string = "{:." + str(num_digits) + "e}"
        return format_string.format(value)
    else:
        format_string = "{:." + str(num_digits) + "f}"
        return format_string.format(value * 100) + "%"


def write_object_to_pickle(file_name, to_pickle_object):
    with open(file_name, 'wb') as handle:
        pickle.dump(to_pickle_object, handle, protocol=pickle.HIGHEST_PROTOCOL)


def get_tissue_name(sample_name):
    splitted = sample_name.split("-")
    if (len(splitted) > 1):
        without_id = splitted[:-1]
        return '-'.join(without_id)
    else:
        return sample_name


def get_tissue_names(data_path):
    tissue_names = []
    for filename in os.listdir(data_path):
        if filename.endswith(".gz"):
            sample_name = filename.split('.')[0]
            tissue_name = get_tissue_name(sample_name)
            tissue_names.append(tissue_name)

    return set(tissue_names)

def get_sample_names(directory):
    sample_names = []
    for filename in os.listdir(directory):
        if filename.endswith(".gz"):
            sample_name = filename.split('.')[0]
            sample_names.append(sample_name)
    return sample_names

def get_num_samples_by_tissue(tissue_name, sample_names):
    count = 0
    for sample_name in sample_names:
        if get_tissue_name(sample_name) == tissue_name:
            count+=1
    return count


def get_tissue_file_dict(groups_file_name="/cs/cbio/jon/groups.csv"):
    df = pd.read_csv(groups_file_name)
    # df.loc[df.group == "Adipocytes", 'include'] = True
    # df = df[df.include == True]
    df = df[~df.group.str.contains(":")]
    grouped_by_tissue = df.groupby(by="group")
    tissue_file_dict = {}
    for name, group in grouped_by_tissue:
        file_list = []
        for index, row in group.iterrows():
            file_list.append(row['name'])
        tissue_file_dict[name] = file_list
    return tissue_file_dict


def get_tissues_file_dict_from_chart(df):#samples_table = pd.read_csv("groups.csv")
    df.loc[df.group == "Adipocytes", 'include'] = True
    df = df[df.include == True]
    df = df[~df.group.str.contains(":")]
    grouped_by_tissue = df.groupby(by="group")
    tissue_file_dict = {}
    for name, group in grouped_by_tissue:
        file_list = []
        for index, row in group.iterrows():
            file_list.append(row['name'])
        tissue_file_dict[name] = file_list
    return tissue_file_dict

def get_count_from_string(letter_to_check, in_string):
    cur_count = 0
    for index, letter in enumerate(in_string):
        if letter == letter_to_check:
            cur_count = "" + in_string[index + 1]
            cur_index = index + 2
            while cur_index < len(in_string):
                cur_char = in_string[cur_index]
                if cur_char.isalpha():
                    break
                cur_index += 1
                cur_count = cur_count + cur_char
            cur_count = int(cur_count)
            break
    return cur_count


def get_unknown_letter_count_from_string(in_string):
    count_dict = {}
    known_letters = {'A', 'C', 'T', 'G'}
    for index, letter in enumerate(in_string):
        if letter.isalpha() and not letter in known_letters:
            cur_count = "" + in_string[index + 1]
            cur_index = index + 2
            while cur_index < len(in_string):
                cur_char = in_string[cur_index]
                if cur_char.isalpha():
                    break
                cur_index += 1
                cur_count = cur_count + cur_char
            cur_count = int(cur_count)
            count_dict[letter] = cur_count
    return count_dict


def is_two_blocks(start_ind, end_ind, beta_arrays):
    min_block_size = 5
    if start_ind > end_ind - (2 * min_block_size):
        return True
    left_avg = 0
    left_count = 0
    global_avg = 0
    global_total = 0
    right_avg = 0
    right_count = 0
    min_variance = 10000
    min_variance_index = -1
    cur_ind = start_ind
    while cur_ind < end_ind:
        for b_array in beta_arrays:
            row = b_array[cur_ind]
            methyl_prop = 0 if row[1] == 0 else row[0] / row[1]
            global_avg += methyl_prop
            global_total += 1
            if (cur_ind < min_block_size):
                left_avg += methyl_prop
                left_count += 1
            else:
                right_avg += methyl_prop
                right_count += 1
        cur_ind += 1

    right_avg = 0 if right_count == 0 else right_avg / right_count
    left_avg = 0 if left_count == 0 else left_avg / left_count
    global_avg = 0 if global_total == 0 else global_avg / global_total

    cur_ind = start_ind + min_block_size
    while cur_ind <= end_ind - min_block_size:
        cur_site_sum = 0
        num_samples = 0
        for b_array in beta_arrays:
            row = b_array[cur_ind]
            methyl_prop = 0 if row[1] == 0 else row[0] / row[1]
            cur_site_sum += methyl_prop
            num_samples += 1
        left_avg = 0 if (left_count + num_samples) == 0 else ((left_avg * left_count) + cur_site_sum) / (
                left_count + num_samples)
        left_count += num_samples
        right_avg = 0 if (right_count - num_samples) == 0 else ((right_avg * right_count) - cur_site_sum) / (
                right_count - num_samples)
        right_count -= num_samples
        inner_ind = start_ind
        diff_left = 0
        diff_right = 0
        diff_global = 0
        diff_count_left = 0
        diff_count_right = 0
        while inner_ind < end_ind:
            cur_site_sum = 0
            num_samples = 0
            for b_array in beta_arrays:
                row = b_array[inner_ind]
                methyl_prop = 0 if row[1] == 0 else row[0] / row[1]
                num_samples += 1
                diff_global += (methyl_prop - global_avg) ** 2
                if (inner_ind <= cur_ind):
                    diff_left += (methyl_prop - left_avg) ** 2
                    diff_count_left += 1
                else:
                    diff_right += (methyl_prop - right_avg) ** 2
                    diff_count_right += 1
            inner_ind += 1

        diff_split = 0 if (diff_count_left + diff_count_right) == 0 else (diff_left + diff_right) / (
                diff_count_left + diff_count_right)
        diff_global = 0 if (diff_count_left + diff_count_right) == 0 else diff_global / (
                diff_count_left + diff_count_right)
        if diff_split < min_variance:
            min_variance = diff_split
            min_variance_index = cur_ind
        cur_ind += 1

    global_delta = (diff_global / 2) + 0.0001

    return diff_global - global_delta > min_variance


def get_is_2_blocks(blocks_df, beta_arrays):
    average_dist_from_half = 0
    total = 0
    avg_dis_from_half_list = []
    for index, row in blocks_df.iterrows():
        start_ind = row.d - 1
        end_ind = row.e - 1
        a_bool = is_two_blocks(start_ind, end_ind, beta_arrays)
        if a_bool:
            average_dist_from_half = 1
        avg_dis_from_half_list.append(a_bool)
    return avg_dis_from_half_list


# input is a homog file line
def line_to_parts(sample_line):
    split_line = sample_line.decode("utf-8")[:-1].split("\t")
    u_sites = int(split_line[5])
    x_sites = int(split_line[6])
    m_sites = int(split_line[7])
    start = int(split_line[1])
    end = int(split_line[2])
    region = split_line[0] + ":" + split_line[1] + "-" + split_line[2]
    return u_sites, x_sites, m_sites, region


def break_to_chunks(start_ind, end_ind, num_chunks):
    step_size = int((end_ind - start_ind) / num_chunks)
    range_list = list(range(start_ind, end_ind, step_size))
    res_list = []
    for idx, val in enumerate(range_list):
        if idx == len(range_list) - 1:
            break
        res_list.append((val, range_list[idx + 1]))
    res_list.append((range_list[-1], end_ind))
    return res_list


class DiscoveryMethod(Enum):
    MIN_UM = 1
    ENTROPY = 2


class UmScore(Enum):
    AVG_KEY = 'min_um_avg'
    DEVIATION_KEY = 'min_um_deviation'


def first_pass_stats_by_key(region_dict, key='min_um'):
    for region, region_map in region_dict.items():
        avg_entropy = region_map[key + '_sum'] / region_map['total']
        region_map[key + '_avg'] = avg_entropy


def second_pass_stats_by_key(region_dict, metric_name='min_um'):
    counter = 0
    for region in region_dict:
        #         if (counter % 100000 == 0):
        #             print("progress")
        counter += 1
        entropy_deviation = region_dict[region][metric_name + '_deviation']
        global_variance_etnropy = entropy_deviation / region_dict[region]['total']
        region_dict[region]['standard_dev_' + metric_name] = global_variance_etnropy ** 0.5


def z_score_by_region_map(global_region_propotion_dict, sampled_region_proportion_dict):
    score_map = {}
    for region, region_map in sampled_region_proportion_dict.items():
        sum_avg = region_map['min_um_avg']
        global_map = global_region_propotion_dict[region]
        global_sum_avg = global_map['min_um_avg']
        global_stnd_dev_sum = global_map['standard_dev_min_um'] + 0.00001
        z_score_sum = 0.0 if sum_avg == 0 else (sum_avg - global_sum_avg) / global_stnd_dev_sum
        score_map[region] = z_score_sum
    return score_map


def save_candidate_file(candidates, file_name="candidates.bed"):
    with open(file_name, "w") as text_file:
        for block in candidates:
            chrome = block.split(':')[0]
            start_ind = block.split(':')[1].split('-')[0]
            end_ind = block.split(':')[1].split('-')[1]
            print(chrome + '\t' + start_ind + '\t' + end_ind, file=text_file)


def save_candidate_bimodal_files_list(candidates, file_name="candidates.bed"):
    with open(file_name, "w") as text_file:
        for block, file_list in candidates:
            chrome = block.split(':')[0]
            start_ind = block.split(':')[1].split('-')[0]
            end_ind = block.split(':')[1].split('-')[1]
            file_list_str = ":".join(file_list)
            print(chrome + '\t' + start_ind + '\t' + end_ind + '\t' + file_list_str, file=text_file)


def save_candidate_tuple_file(candidates, file_name="candidates.bed"):
    with open(file_name, "w") as text_file:
        for block, score in candidates:
            chrome = block.split(':')[0]
            start_ind = block.split(':')[1].split('-')[0]
            end_ind = block.split(':')[1].split('-')[1]
            print(chrome + '\t' + start_ind + '\t' + end_ind + '\t' + str(score), file=text_file)


def awk_program(methylated_reads):
    if methylated_reads:
        var_name = "methyl_prop"
    else:
        var_name = "unmethyl_prop"

    return """
        awk '{
          for (i=1;i<=NF;i++)
              {
                if ($i ~ /YI/) {
                    split($i,methyl_info,\":\");
                    split(methyl_info[3],methyl_counts,\",\");
                    total_cpgs = methyl_counts[1] + methyl_counts[2];
                    if (total_cpgs > 5){
                        unmethyl_prop = methyl_counts[2] / total_cpgs;
                        methyl_prop = methyl_counts[1] / total_cpgs;
                        if (%s > 0.75){
                            print $0;
                        }
                    }
                }
              }
         }'
     """ % (var_name)


def create_pileup_file(candidates_file, bam_file):
    return 0


def get_pileup_file(filename, file_sufix=".vcf"):
    filename_no_extension = os.path.splitext(filename)[0]
    out_file = filename_no_extension + "_pileup" + file_sufix
    return out_file


def get_tissue_specific_bimodal_blocks_dict_name(tissue_name):
    return tissue_name + "_bimodal_blocks_dict.pkl"


def get_tissue_specific_dict_name(tissue_name):
    return tissue_name + "_region_sum_dict.pkl"


def run_command(cmd, stderr_out=subprocess.STDOUT):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=stderr_out)
    output, error = process.communicate()
    return output.decode("utf-8")


class SNPResult:
    def __init__(self, glbl_letters, glbl_counts, m_letters, m_counts, u_letters, u_counts, pos, samplename):
        self.glbl_letters = glbl_letters
        self.glbl_counts = glbl_counts
        self.m_letters = m_letters
        self.m_counts = m_counts
        self.u_letters = u_letters
        self.u_counts = u_counts
        self.pos = pos
        self.samplename = samplename


def find_gene_names(block, reference_path="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/",
                    annotation_file="annotations.bed.gz", block_window=100):
    new_region = add_space_to_block(block, block_window)
    find_gene_names_cmd = ["tabix", reference_path + annotation_file, new_region]
    result = subprocess.run(find_gene_names_cmd, stdout=subprocess.PIPE)
    res = result.stdout.decode('utf-8').split("\n")
    gene_set = []
    for line in res:
        if len(line) > 0:
            gene_name = line.split("\t")[-1]
            gene_set.append(gene_name)
    return set(gene_set)


def get_gene_lines(block, reference_path="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/",
                    annotation_file="annotations.bed.gz", block_window=100):
    new_region = add_space_to_block(block, block_window)
    find_gene_names_cmd = ["tabix", reference_path + annotation_file, new_region]
    result = subprocess.run(find_gene_names_cmd, stdout=subprocess.PIPE)
    res = result.stdout.decode('utf-8').split("\n")
    gene_set = []
    for line in res:
        if len(line) > 0:
            tokens = line.split("\t")
            gene_name = tokens[-1]
            start = int(tokens[1])
            end = int(tokens[2])
            gene_set.append((gene_name, start, end))
    return set(gene_set)


def find_lncRNA_names(blocks, reference_path="/cs/cbio/jon/lncRNA/",
                    annotation_file="hg19.lncRNAs.sorted.bed.gz", block_window=1000):
    found_lnr_rnas = {}
    for block in blocks:
        new_region = add_space_to_block(block, block_window)
        find_gene_names_cmd = ["tabix", reference_path + annotation_file, new_region]
        result = subprocess.run(find_gene_names_cmd, stdout=subprocess.PIPE)
        res = result.stdout.decode('utf-8').split("\n")
        for line in res:
            if len(line) > 0:
                tokens = line.split("\t")
                gene_name = tokens[3]
                cur_start = int(tokens[1])
                cur_end = int(tokens[2])
                cur_chrom = tokens[0]
                if cur_chrom in found_lnr_rnas:
                    cur_lnc_rna_dict = found_lnr_rnas[cur_chrom]
                else:
                    cur_lnc_rna_dict = {}
                if gene_name in cur_lnc_rna_dict:
                    gene_start_end = cur_lnc_rna_dict[gene_name]
                    gene_start = gene_start_end[0]
                    gene_end = gene_start_end[1]
                    gene_start = min(gene_start, cur_start)
                    gene_end = max(gene_end, cur_end)
                    cur_lnc_rna_dict[gene_name] = (gene_start, gene_end)
                else:
                    cur_lnc_rna_dict[gene_name] = (cur_start, cur_end)
                found_lnr_rnas[cur_chrom] = cur_lnc_rna_dict
    return found_lnr_rnas


def get_chroms():
    cpg_path = '/cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.bed.gz'
    bp_cpg_mat = pd.read_csv(cpg_path, sep="\t", header=None, names=["chr", "start", "cpg_start"])
    chroms = list(set(bp_cpg_mat.chr))
    return chroms


def find_gene_start_end(gene_name,
                    annotation_file="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/annotations.bed.gz"):
    find_gene_cmd = " ".join(["zcat", annotation_file, "| grep", gene_name])

    res = run_command(find_gene_cmd)
    res = res.split("\n")
    if len(res) > 0 and res[0] != '':
        start = 0
        end = 0
        chrom = ""
        for line in res:
            if len(line) > 0:
                # line = line.decode()
                tokens = line.split("\t")
                cur_gene_name = tokens[-1]
                if cur_gene_name == gene_name:
                    chrom = tokens[0]
                    cur_start = int(tokens[1])
                    cur_end = int(tokens[2])
                    if cur_start < start or start == 0:
                        start = cur_start
                    if cur_end > end:
                        end = cur_end
        return chrom, start, end
    else:
        return None, None, None


def genes_in_window_with_distance(block, gene_name=None, reference_path="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/",
                    annotation_file="annotations.bed.gz"):
    find_gene_names_cmd = ["tabix", reference_path + annotation_file, block]

    result = subprocess.run(find_gene_names_cmd, stdout=subprocess.PIPE)
    res = result.stdout.decode('utf-8').split("\n")
    found_genes = []
    seen_genes = set()
    for line in res:
        if len(line) > 0:
            tokens = line.split("\t")
            cur_gene_name = tokens[-1]
            cur_anno_type = tokens[3]
            if gene_name != cur_gene_name and cur_anno_type == "exon":
                if cur_gene_name not in seen_genes:
                    cur_start = int(tokens[1])
                    cur_end = int(tokens[2])
                    found_genes.append((cur_gene_name, cur_start, cur_end))
                    seen_genes.add(cur_gene_name)
    return found_genes


def find_nearest_genes(gene_name, reference_path="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/",
                    annotation_file="annotations.bed.gz", num_neighbors=3):
    block_window = 1e5
    res = find_gene_start_end(gene_name)
    if res is not None:
        chrom = res[0]
        start = res[1]
        end = res[2]
        if chrom == "":
            return []
        block = "{}:{}-{}".format(chrom, str(start), str(end))
        new_region = add_space_to_block(block, block_window)
        found_near_genes = []
        buffer = 0
        while len(found_near_genes) < num_neighbors:
            found_near_genes = genes_in_window_with_distance(new_region, gene_name)
            buffer = buffer + 50000
            new_region = add_space_to_block(block, block_window + buffer)
        genes_distances = []
        for cur_gene, cur_start, cur_end in found_near_genes:
            distance = min(abs(cur_end - start), abs(cur_end - end), abs(cur_start - start), abs(cur_start - end))
            genes_distances.append((cur_gene, distance))
        genes_distances = list(sorted(genes_distances, key=lambda item: item[1]))
        genes_distances = genes_distances[0:num_neighbors]
        return [gn for gn, _ in genes_distances]
    else:
        return []


def find_gene_block(gene_name, reference_path="/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/references/hg19/",
                    annotation_file="annotations.bed.gz"):
    annotation_gene_lookup_cmd = "zcat {}{}| grep {}".format(reference_path, annotation_file, gene_name)
    result = run_command(annotation_gene_lookup_cmd)
    res = result.split("\n")
    start_site = 51244562 #last possible start_site
    end_site = 0
    chrom = ""
    for lin in res:
        tokens = lin.split("\t")
        cur_gene_name = tokens[-1]
        if cur_gene_name == gene_name:
            chrom = tokens[0]
            cur_start_site = int(tokens[1])
            cur_end_site = int(tokens[2])
            if cur_start_site < start_site:
                start_site = cur_start_site
            if cur_end_site > end_site:
                end_site = cur_end_site
    return chrom, start_site, end_site


def add_space_to_block(block, block_window=100):
    region_inds = block.split(":")[1].split("-")
    new_region = block.split(":")[0] + ":" + str(int(int(region_inds[0]) - block_window)) + "-" + str(
        int(int(region_inds[1]) + block_window))
    return new_region