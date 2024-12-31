import optparse
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

from asm_utils import break_to_chunks, run_command
from subprocess import Popen, PIPE
from multiprocessing import Pool
from scipy.stats import chisquare, distributions
import numpy as np

from qvalue import q_estimate

coverage_threshold = 5
bimodal_proportion = 0.27
homozygous_threshold = 0.8


class SNP:
    hetero_m_side_ref_count = 0
    hetero_u_side_ref_count = 0
    hetero_m_side_alt_count = 0
    hetero_u_side_alt_count = 0
    m_block_uniform_m_side_homo_ref_count = 0
    m_block_uniform_u_side_homo_ref_count = 0
    m_block_uniform_m_side_homo_alt_count = 0
    m_block_uniform_u_side_homo_alt_count = 0
    m_block_u_side_uniform_hetero_count = 0
    m_block_m_side_uniform_hetero_count = 0
    u_block_uniform_m_side_homo_ref_count = 0
    u_block_uniform_u_side_homo_ref_count = 0
    u_block_uniform_m_side_homo_alt_count = 0
    u_block_uniform_u_side_homo_alt_count = 0
    u_block_u_side_uniform_hetero_count = 0
    u_block_m_side_uniform_hetero_count = 0
    b_block_uniform_m_side_homo_ref_count = 0
    b_block_uniform_u_side_homo_ref_count = 0
    b_block_M_hetero = 0
    b_block_U_hetero = 0
    M_hetero_letters = {}
    U_hetero_letters = {}
    b_block_uniform_ref_count = 0
    b_block_uniform_alt_count = 0
    count_danger_side = 0
    count_safe_side = 0
    ref_letter = ''
    alt_letter = ''
    files_set = set()

    def __init__(self):
        self.files_set = set()
        self.M_hetero_letters = {}
        self.U_hetero_letters = {}

    def add_to_files_set(self, filename):
        self.files_set.add(filename)


def get_u_m_prop(u, x, m):
    total = u + x + m
    u_prop = u / total if total > 0 else 0
    m_prop = m / total if total > 0 else 0
    return u_prop, m_prop


def is_bimodal(u_prop, m_prop):
    return min(u_prop, m_prop) > bimodal_proportion


def is_hetero(letters_array, counts_array):
    # hetero_threshold = 4
    total = sum(counts_array)
    prop_1 = counts_array[0] / total if total > 0 else 0
    prop_2 = counts_array[1] / total if total > 0 else 0
    return letters_array[0] != letters_array[1] and (prop_1 >= 0.2 and prop_2 >= 0.2 and (counts_array[0] + counts_array[1]) > 3)
    # return letters_array[0] != letters_array[1] and counts_array[0] > hetero_threshold and counts_array[1] > hetero_threshold


def find_meth_letters_and_counts(snp_info_array, meth_status):
    for ind, el in enumerate(snp_info_array):
        if el == meth_status:
            return snp_info_array[ind + 1], snp_info_array[ind + 2]
    return "./.", "0/0"


def get_main_letter(letters, counts):
    letters = letters.split("/")
    counts = [int(el) for el in counts.split("/")]
    total = counts[0] + counts[1]
    prop_1 = counts[0] / total if total > 0 else 0
    prop_2 = counts[1] / total if total > 0 else 0
    if total < 2:
        return ".", 0, total
    if prop_1 > homozygous_threshold:#before >=
        return letters[0], total, total
    elif prop_2 > homozygous_threshold:#before >=
        return letters[1], total, total
    elif (counts[0] > 3 and counts[1] > 3) or (prop_1 >= 0.2 and prop_2 >= 0.2 and (counts[0] + counts[1]) > 3):
        return "H", total, total
    else:
        return ".", 0, total


def get_snps_info_extraction(snp_el):
    snp_info_array = snp_el.split(",")
    cur_snp_pos_bp = snp_info_array[0]
    glbl_letters = snp_info_array[1]
    glbl_counts = snp_info_array[2]
    ref_letter = glbl_letters.split("/")[0]
    alt_letter = glbl_letters.split("/")[1]
    main_letter, total_count, _ = get_main_letter(glbl_letters, glbl_counts)
    return snp_info_array, cur_snp_pos_bp, glbl_letters, glbl_counts, ref_letter, alt_letter, main_letter, total_count


def get_u_m_split_checks(snp_info_array):
    m_letters, m_counts = find_meth_letters_and_counts(snp_info_array, "M")
    u_letters, u_counts = find_meth_letters_and_counts(snp_info_array, "U")
    m_letter, m_total, m_real_total = get_main_letter(m_letters, m_counts)
    u_letter, u_total, u_real_total = get_main_letter(u_letters, u_counts)
    return m_letter, m_total, m_real_total, u_letter, u_total, u_real_total


def get_u_m_counts(snp_info_array):
    m_letters, m_counts = find_meth_letters_and_counts(snp_info_array, "M")
    u_letters, u_counts = find_meth_letters_and_counts(snp_info_array, "U")
    return m_letters, m_counts, u_letters, u_counts


def choose_blocks_by_fdr_bh(pvals, blocks, alpha=0.1, return_only_selected=True):
    rejected_list, corrected_p_vals, _, _ = multipletests(pvals, alpha=alpha, method='fdr_bh')
    index = 0
    if return_only_selected:
        for rejected in rejected_list:
            if not rejected:
                break
            index += 1
        if index > 0:
            return corrected_p_vals[0:index], blocks[0:index], pvals[0:index]
        else:
            return [], [], []
    else:
        return corrected_p_vals, blocks, pvals


def choose_blocks_for_cell_heterogeneiety(pvals, blocks, threshold=0.5):
    index = 0
    for p_val in pvals:
        if p_val > threshold:
            break
        index += 1
    if index > 0:
        return pvals[index:], blocks[index:]
    else:
        return [], []


def choose_blocks_by_q_value(pvals, blocks, significance=0.05):
    q_values = q_estimate(np.array(pvals))
    index = 0
    for q_val in q_values:
        if q_val > significance:
            break
        index += 1
    if index > 0:
        return pvals[0:index], blocks[0:index]
    else:
        return [], []


def calculate_chi_squared(m_counts, u_counts, glbl_counts, test_name="fisher", psuedo_count=0):
    counts_array = np.zeros([2, 2])
    counts_array[0, 0] = m_counts[0]
    counts_array[0, 1] = m_counts[1]
    counts_array[1, 0] = u_counts[0]
    counts_array[1, 1] = u_counts[1]
    counts_array = counts_array + psuedo_count
    
    if test_name == "chi_square":
        b_probs = np.zeros(counts_array.shape[1])
        
        total = counts_array[0, 0] + counts_array[0, 1] + counts_array[1, 0] + counts_array[1, 1]
        for i in range(counts_array.shape[1]):
            b_probs[i] = sum(counts_array[:, i]) / total
        # b_probs[0] = glbl_counts[0] / sum(glbl_counts)
        # b_probs[1] = 1 - b_probs[0]
    
        expected_array = np.zeros([2, 2])
        for i in range(2):
            for j in range(2):
                expected_array[i, j] = b_probs[j] * sum(counts_array[i, :])
        if b_probs[0] == 1 or b_probs[1] == 1:
            return 0, 1
        # counts_array = counts_array + 1
        # expected_array = expected_array + 1
        chi_squared, p_val = chisquare(counts_array.flatten(), f_exp=expected_array.flatten())
        return chi_squared, p_val
    elif test_name == "boschloo":
        res = stats.boschloo_exact(counts_array, alternative="two-sided")
        return res.statistic, res.pvalue
    elif test_name == "fisher":
        oddsr, p = stats.fisher_exact(counts_array, alternative="two-sided")
        return oddsr, p
    else:
        raise NotImplementedError


class SnpProcessor:
    seent_it = False
    snp_dict = {}
    sd_asm_bimodal_list = []
    sd_asm_list = []
    perfect_sd_list = []
    imprinting_list = []
    imprinting_res_list = []
                            
    def choose_snp_with_most_coverage(self, snps_list):
        max_snp = None
        max_total = 0
        for snp_el in snps_list:
            snp_info_array, cur_snp_pos_bp, glbl_letters, glbl_counts, ref_letter, alt_letter, letter, total = \
                get_snps_info_extraction(snp_el)
            # if cur_snp_pos_bp == "9826962":
            #     x = 0
            glbl_letters = glbl_letters.split("/")
            glbl_counts = [int(x) for x in glbl_counts.split("/")]
            if is_hetero(glbl_letters, glbl_counts):
                if total > max_total:
                    max_total = total
                    max_snp = snp_el
        return max_snp


    def iterate_over_snps_list_chi_Squared_sm(self, snps_list):
        snps_tokens = snps_list.split(";")
        for snp_el in snps_tokens:
            snp_info_array, cur_snp_pos_bp, glbl_letters, glbl_counts, ref_letter, alt_letter, letter, total = \
                get_snps_info_extraction(snp_el)
            # snp_key = chrom + ":" + str(cur_snp_pos_bp)
            chi_squared_sum = 0
            deg_of_freedom = 0
            if is_hetero(glbl_letters, glbl_counts):
                if total > coverage_threshold:
                    m_letters, m_counts, u_letters, u_counts = get_u_m_counts(snp_info_array)
                    m_letters = m_letters.split("/")
                    m_counts = [int(x) for x in m_counts.split("/")]
                    u_letters = u_letters.split("/")
                    u_counts = [int(x) for x in u_counts.split("/")]
                    if m_letters == u_letters and is_hetero(m_letters, m_counts) and is_hetero(u_letters, u_counts):
                        chi_squared, p_val = calculate_chi_squared(m_counts, u_counts)
                        chi_squared_sum += chi_squared
                        deg_of_freedom += 1
            p = distributions.chi2.sf(chi_squared_sum, deg_of_freedom)

    def iterate_over_snps_get_p_val_hetero_snp(self, snps_list):
        snp_el = self.choose_snp_with_most_coverage(snps_list)
        if snp_el:
            snp_info_array, snp_pos, _, glbl_counts, _, _, _, _ = get_snps_info_extraction(snp_el)
            # if snp_pos == "119233110":
            #     x = 0
            m_letters, m_counts, u_letters, u_counts = get_u_m_counts(snp_info_array)
            m_letters = m_letters.split("/")
            m_counts = [int(x) for x in m_counts.split("/")]
            u_letters = u_letters.split("/")
            u_counts = [int(x) for x in u_counts.split("/")]
            glbl_counts = [int(x) for x in glbl_counts.split("/")]
            if m_letters == u_letters and not is_hetero(m_letters, m_counts) and not is_hetero(u_letters, u_counts) and np.sum(m_counts) > 3 and np.sum(u_counts) > 3:
                chi_squared, p_val = calculate_chi_squared(m_counts, u_counts, glbl_counts)
                return chi_squared, p_val, snp_el
        return None, None, None


    def iterate_over_blocks(self, file_name, block, in_memory=True): # in_memory=False should not be used with many threads . it crashes
        file_iterators_list = []
        # chr4:69978091-69978091
        # block = "chr19:10706400-10708400"
        tabix_cmd = ["tabix", file_name, block]
        p = Popen(tabix_cmd, stdout=PIPE, bufsize=1,
                   universal_newlines=True)
        # f = gz.open(directory + filename, 'r')
        # for f in file_iterators_list:
        #     cur_line = f.readline()
        should_continue = True
        p_val_list = []
        while should_continue:
            cur_line = p.stdout.readline().rstrip()#.decode().rstrip()
            if len(cur_line) > 3:
                cur_tokens = cur_line.split("\t")
                if len(cur_tokens) > 10:
                    snps_list = cur_tokens[-1].split(";")
                    chrom = cur_tokens[0]
                    region = chrom + ":" + cur_tokens[1] + "-" + cur_tokens[2]
                    u = int(cur_tokens[5])
                    x = int(cur_tokens[6])
                    m = int(cur_tokens[7])
                    u_prop, m_prop = get_u_m_prop(u, x, m)
                    filename = p.args[1].split("/")[-1].split(".")[0]
                    # TODO handle not counting snps multiple times in the same file from different blocks
                    if is_bimodal(u_prop, m_prop):
                        chi_squared, p_val, snp_el = self.iterate_over_snps_get_p_val_hetero_snp(snps_list)
                        if p_val is not None:
                            p_val_list.append((p_val, (cur_tokens[:-1], snp_el)))
            else:
                should_continue = False
                break
        p.wait()
            # f.close()
        p_val_list = sorted(p_val_list, key=lambda item: item[0])

        return p_val_list


def iterate_over_blocks_wrapper(file_name, block, block_override=None):
    print(block)
    snp_processor = SnpProcessor()
    if block_override:
        block = block_override
    return snp_processor.iterate_over_blocks(file_name, block)


def run_process(file_name, cell_heterogeniety, num_threads=16, block=None, block_file="/cs/zbio/jrosensk/block_files_2/homog/all_blocks/all_blocks.sorted.tsv.gz"):
    if block:
        chrom = block.split(":")[0]
        chrom_list = [chrom]
    else:
        chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                  "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                  "chr22", "chrX", "chrY"]
        # chrom_list = ["chr4"]

    params_list = []

    if block:
        params_list = [(file_name, None, block)]
    else:
        for chrom in chrom_list:
            tabix_command = "tabix {} {}".format(block_file, chrom)
            block_output = run_command(tabix_command)
            block_lines = block_output.split("\n")
            if len(block_lines[-1]) <= 1:
                block_lines = block_lines[:-1]
            region_start = int(block_lines[0].split("\t")[1])
            region_end = int(block_lines[-1].split("\t")[2])
            chunks_list = break_to_chunks(region_start, region_end, 3)
            for chunk in chunks_list:
                region = chrom + ":" + str(chunk[0]) + "-" + str(chunk[1])
                params_list.append((file_name, region, None))
    # params_list.append(("/cs/zbio/jrosensk/homog_files/", "chr4:140656829-140656859"))

    with Pool(num_threads) as p:
        res = p.starmap(iterate_over_blocks_wrapper, params_list)
        p.close()
        p.join()
    p_val_list = [item for sublist in res for item in sublist]
    p_val_list = sorted(p_val_list, key=lambda item: item[0])
    res_list = [list(t) for t in zip(*p_val_list)]
    sorted_p_vals = res_list[0]
    blocks_and_snp = res_list[1]
    if cell_heterogeniety:
        p_vals, res_blocks = choose_blocks_for_cell_heterogeneiety(sorted_p_vals, blocks_and_snp, 0.8)
    else:
        p_vals, res_blocks, _ = choose_blocks_by_fdr_bh(sorted_p_vals, blocks_and_snp)
    return p_vals, res_blocks


if __name__ == '__main__':
    # file_name_new = "/cs/zbio/jrosensk/efrat_derived/homog/homog_aligned/Aorta-Endothel-Z0000043G.snps.gz"
    # file_name = "/cs/zbio/jrosensk/homog_files/Bladder-Epithelial-Z0000043F.snps.gz"
    parser = optparse.OptionParser()
    parser.add_option('--in_file',#
                      default="/cs/zbio/jrosensk/block_files_smoothed/homog/homog_aligned/Blood-Granulocytes-Z000000UD.snps.gz")
    parser.add_option('--block_file',#
                      default="/cs/zbio/jrosensk/road_map_block_files/homog/all_blocks/all_blocks.sorted.tsv.gz")
    parser.add_option('--out_file',
                      default="/cs/zbio/jrosensk/road_map_block_files/homog/homog_aligned/out.asm")
    parser.add_option('--cell_heterogeneity', action="store_true")
    parser.add_option('--num_threads',
                      default=16)
    options, arguments = parser.parse_args()
    file_name = options.in_file
    if options.out_file is None:
        out_file = file_name.split(".")[0] + ".asm"
    else:
        out_file = options.out_file
    block = None#"chr1:3339883-3340235"
    # print(file_name)
    p_vals, res_blocks_and_snps = run_process(file_name, options.cell_heterogeneity, num_threads=int(options.num_threads),
                                              block=block, block_file=options.block_file)
    with open(out_file, 'w') as f_out:
        for p_val, (block, snp_el) in zip(p_vals, res_blocks_and_snps):
            # chrom = block.split(":")[0]
            # range = block.split(":")[-1]
            # range_tokens = range.split("-")
            # start = range_tokens[0]
            # end = range_tokens[1]
            snp_info_array, cur_snp_pos_bp, glbl_letters, glbl_counts, ref_letter, alt_letter, letter, total = \
                get_snps_info_extraction(snp_el)
            m_letters, m_counts, u_letters, u_counts = get_u_m_counts(snp_info_array)
            block_line = "\t".join(block)
            f_out.write(f"{block_line}\t{p_val:,.3e}\t{cur_snp_pos_bp}\t{glbl_letters}\t{glbl_counts}\t{m_counts}\t{u_counts}\n")
