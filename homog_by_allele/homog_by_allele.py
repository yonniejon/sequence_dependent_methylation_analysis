#!/usr/bin/python3 -u

"""
Allele-Specific Methylation Homogeneity Analyzer
------------------------------------------------
This script quantifies DNA methylation patterns across specific alleles by
integrating WGBS data (BAM files) with known SNP locations. It calculates the
distribution of unmethylated, partially methylated, and fully methylated
fragments (UXM) for each allele, enabling the detection of Allele-Specific
Methylation (ASM) influenced by genetic variants.

Key Capabilities:
- Parallel Allele Processing: Distributes genomic regions across multiple
  threads to analyze BAM files against SNP libraries efficiently.
- SNP-Guided Analysis: Uses tabix-indexed SNP files to distinguish between
  reference and alternative alleles within read fragments.
- Flexible Filtering: Supports custom mapping quality (MAPQ) thresholds,
  SAM flag exclusions, and minimum CpG site requirements per read.
- Fragment Classification: Categorizes reads into "U" (unmethylated), "X"
  (intermediate), or "M" (methylated) bins based on user-defined thresholds.
"""

import os
import os.path as op
import subprocess
import shlex
import datetime
from multiprocessing import Pool
import argparse
from init_genome import chromosome_order
from utils_wgbs import IllegalArgumentError, match_maker_tool, eprint, validate_local_exe
from bam2pat import add_args, subprocess_wrap, validate_bam, is_pair_end, MAPQ, extend_region
from genomic_region import GenomicRegion


FILE_SUF = 'allele.homog'

# Minimal Mapping Quality to consider.
# 10 means include only reads w.p. >= 0.9 to be mapped correctly.
# And missing values (255)

homog_by_allele_tool = 'cpp/SnpMethHist'


def proc_chr(input_path, out_path_name, region, genome, paired_end, ex_flags, mapq, min_cpg,
             snps_file, regions_file, in_flags, rate_cmd):
    """ Convert a temp single chromosome file, extracted from a bam file,
        into a sam formatted (no header) output file."""

    # Run patter tool 'bam' mode on a single chromosome

    out_path = out_path_name + f'.output.{FILE_SUF}'
    bed_file = None

    # use samtools to extract only the reads from 'chrom'
    # flag = '-f 3' if paired_end else ''
    if in_flags is None:
        in_flags = '-f 3' if paired_end else ''
    else:
        in_flags = f'-f {in_flags}'
    cmd = "samtools view {} {} -q {} -F {} {} -M -L {}".format(input_path, region, mapq, ex_flags, in_flags, regions_file)
    if bed_file is not None:
        cmd += f" -M -L {bed_file} | "
    else:
        cmd += "| "
    if paired_end:
        # change reads order, s.t paired reads will appear in adjacent lines
        cmd += f'{match_maker_tool} | '
    cmd += f'{homog_by_allele_tool} {genome.dict_path} {snps_file} /cs/cbio/jon/projects/SnpCount/blacklist.bed.gz'
    if min_cpg is not None:
        cmd += f' --min_cpg {str(min_cpg)}'
    if rate_cmd:
        cmd += rate_cmd

    cmd += f'> {out_path}'

    # print(cmd)
    subprocess_wrap(cmd, False)
    return out_path


class HomogByAllele:
    def __init__(self, args, bam_path):
        self.args = args
        self.out_dir = args.out_dir
        self.bam_path = bam_path
        self.debug = args.debug
        self.gr = GenomicRegion(args, genome_name=args.genome)
        self.validate_input()

    def validate_input(self):

        # validate bam path:
        validate_bam(self.bam_path)

        # validate output dir:
        if not (op.isdir(self.out_dir)):
            raise IllegalArgumentError('Invalid output dir: {}'.format(self.out_dir))

    def set_regions(self):
        # if user specified a region, just use it
        if self.gr.region_str:
            return [self.gr.region_str]

        bam_chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']# output.decode()[:-1].split('\n')

        # get all chromosomes from the reference genome:
        ref_chroms = self.gr.genome.get_chroms()
        # intersect the chromosomes from the bam and from the reference
        intersected_chroms = list(set(bam_chroms) & set(ref_chroms))

        if not intersected_chroms:
            msg = '[wt acc] Failed retrieving valid chromosome names. '
            msg += 'Perhaps you are using a wrong genome reference. '
            msg += 'Try running:\n\t\twgbstools set_default_ref -ls'
            msg += '\nMake sure the chromosomes in the bam header exists in the reference fasta'
            eprint(msg)
            raise IllegalArgumentError('Failed')

        return list(sorted(intersected_chroms, key=chromosome_order))  # todo use the same order as in ref_chroms instead of resorting it

    def intermediate_bam_file_view(self, name):
        return '<(samtools view {})'.format(name)

    def process_substitute(self, cmd):
        return '<({})'.format(cmd)

    def start_threads(self):
        """ Parse each chromosome file in a different process,
            and concatenate outputs to pat and unq files """
        print(datetime.datetime.now().isoformat() + ": *** starting processing of each chromosome")
        base_name_no_ext = ".".join(op.basename(self.bam_path).split(".")[:-1])
        name = op.join(self.out_dir, base_name_no_ext)
        rate_cmd = ' -r '
        if self.args.thresholds:
            rate_cmd += f'0,{self.args.thresholds},1'
        else:
            raise IllegalArgumentError("Must provide argument 'thresholds'")
        if self.gr.region_str is None:
            final_path = name + f".{FILE_SUF}"
            params = []
            for c in self.set_regions():
                out_path_name = name + '_' + c
                params.append((self.bam_path, out_path_name, c, self.gr.genome,
                        is_pair_end(self.bam_path, self.gr.genome), self.args.exclude_flags,
                        self.args.mapq, self.args.min_cpg, self.args.snps_file, self.args.regions_file, self.args.include_flags, rate_cmd))
            p = Pool(self.args.threads)
            res = p.starmap(proc_chr, params)
            p.close()
            p.join()
        else:
            region_str_for_name = self.gr.region_str.replace(":", "_").replace("-", "_")
            final_path = name + f".{region_str_for_name}" + f".{FILE_SUF}"
            out_path_name = name + '_' + "1"
            res = [proc_chr(self.bam_path, out_path_name, self.gr.region_str, self.gr.genome,
                            is_pair_end(self.bam_path), self.args.exclude_flags, self.args.mapq,
                            self.args.min_cpg, self.args.snps_file, self.args.regions_file, self.args.include_flags, rate_cmd)]
        print('Completed all chromosomes')
        if None in res:
            print('threads failed')
            return

        print(datetime.datetime.now().isoformat() + ": finished processing each chromosome")
        # Concatenate chromosome files

        cmd = f"cat " + ' '.join([p for p in res]) + f" | sort -k1,1 -k2,2n -o {final_path} && bgzip {final_path} && tabix -Cf -b 2 -e 2 {final_path}.gz"
        print(datetime.datetime.now().isoformat() + ': starting cat of files')
        subprocess_wrap(cmd, False)
        print(datetime.datetime.now().isoformat() + ": finished cat of files")

        list(map(os.remove, [l for l in res]))


def add_cpg_args(parser):
    parser.add_argument('--snps_file',  default=None,
                        help='A (tabix indexed) tab seperated file containing SNP locations and ref/alt letters. One row'
                             'example is:\nchr11\t2021164\tG\tT.')
    parser.add_argument('--regions_file', default=None,
                        help='region file to filter bam.')
    parser.add_argument('--thresholds', '-t',
                        help='UXM thresholds, LOW,HIGH. E.g, "0.3334,0.666".\n')
    return parser


def main():
    help_summary = """
        Allele-Specific Methylation (ASM) Analyzer

        This tool integrates WGBS BAM files with SNP data to quantify methylation 
        patterns (UXM) on a per-allele basis. It determines if unmethylated (U), 
        intermediate (X), or methylated (M) fragments are linked to specific 
        genetic variants.

        The pipeline:
        1. Extracts reads from BAM files using samtools.
        2. Filters by MAPQ, SAM flags, and SNP-containing regions.
        3. Categorizes reads into U(nmethylated) (mi)X(ed) or M(ethylated), UXM, bins based on provided thresholds.
        4. Outputs a tabix-indexed file containing allele-specific counts.
        """
    parser = argparse.ArgumentParser(
        description=help_summary,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example: python3 homog_by_allele.py --snps_file snps.txt.gz -t 0.35,0.65 sample.bam"
    )
    parser = add_args(parser)
    parser = add_cpg_args(parser)
    args = parser.parse_args()
    validate_local_exe(homog_by_allele_tool)
    for bam in args.bam:
        if not validate_bam(bam):
            eprint(f'[wt add_cpg_counts] Skipping {bam}')
            continue
        HomogByAllele(args, bam).start_threads()


if __name__ == '__main__':
    main()
