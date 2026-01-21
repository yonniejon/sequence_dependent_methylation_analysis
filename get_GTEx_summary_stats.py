import pandas as pd
import pysam
import argparse
import sys
import os
import time
import tempfile

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
# EBI eQTL Catalogue Metadata File (provided locally)
LOCAL_METADATA_FILE = "helper_data/tabix_ftp_paths.tsv"

# Mapping simplified tissue names to EBI Tissue IDs (sample_group)
TISSUE_MAP = {
    # Custom Mappings (User -> Valid GTEx 'sample_group')
    "Colon-Ep": "colon_sigmoid",
    "Heart-Cardio": "heart_left_ventricle",
    "Kidney-Ep": "kidney_cortex",
    "Liver-Hep": "liver",
    "Pancreas-Acinar": "pancreas",
    "Prostate-Ep": "prostate",
    "Small-Int-Ep": "small_intestine",
    "Gastric-Ep": "stomach",
    "Adipocytes": "adipose_subcutaneous",
    "Blood-Mono+Macro": "blood",
    "Thyroid-Ep": "thyroid",
    "Lung-Ep-Alveo": "lung",
    "Skeletal-Musc": "muscle",
    "Endothel": "artery_aorta",
    "Neuron": "brain_cerebellum",
    "Breast-Basal-Ep": "breast",
    "Esophagus-Ep": "esophagus_gej",
    "Epid-Kerat": "skin_sun_exposed",

    # Standard GTEx Names (Mapped to snake_case sample_group)
    "Adipose_Subcutaneous": "adipose_subcutaneous",
    "Whole_Blood": "blood",
    "Thyroid": "thyroid",
    "Lung": "lung",
    "Muscle_Skeletal": "muscle",
    "Skin_Sun_Exposed_Lower_leg": "skin_sun_exposed"
}

# Hardcoded Column Indices based on EBI eQTL Catalogue Standard Format
# Sample: ENSG... 20 46120612 T C ... 318 0.278 0.793 -0.010 0.039 SNP 967 1340 ...
COL_IDX = {
    'molecular_trait_id': 0,  # Gene ID
    'chromosome': 1,
    'position': 2,
    'maf': 7,
    'pvalue': 8,
    'beta': 9,
    'se': 10,
    'an': 13  # Allele Number (Total alleles, ~ 2 * N)
}


class CaptureStdErr:
    """
    Context manager to capture C-level stderr output (fd 2),
    which is where htslib prints errors that don't raise Python exceptions.
    """

    def __enter__(self):
        sys.stderr.flush()
        self.err_file = tempfile.TemporaryFile(mode='w+')
        self.old_stderr_fd = os.dup(sys.stderr.fileno())
        os.dup2(self.err_file.fileno(), sys.stderr.fileno())
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.flush()
        os.dup2(self.old_stderr_fd, sys.stderr.fileno())
        os.close(self.old_stderr_fd)
        self.err_file.seek(0)
        self.output = self.err_file.read()
        self.err_file.close()


def get_dataset_url(tissue_id):
    """
    Reads the local metadata TSV to find the HTTP path for the summary stats file.
    """
    print(f"Resolving file URL for tissue: {tissue_id}...")

    if not os.path.exists(LOCAL_METADATA_FILE):
        print(f"[Error] Metadata file '{LOCAL_METADATA_FILE}' not found.")
        sys.exit(1)

    try:
        df = pd.read_csv(LOCAL_METADATA_FILE, sep='\t')
        subset = df[
            (df['study_label'] == 'GTEx') &
            (df['sample_group'] == tissue_id) &
            (df['quant_method'] == 'ge')
            ]

        if subset.empty:
            print(f"[Error] No 'ge' dataset found for '{tissue_id}' in the local metadata.")
            sys.exit(1)

        ftp_url = subset.iloc[0]['ftp_path']
        http_url = ftp_url.replace("ftp://", "http://")
        return http_url

    except Exception as e:
        print(f"[Error] Failed to read metadata file: {e}")
        sys.exit(1)


def get_remote_tabix_handle(url):
    """
    Returns a pysam.TabixFile handle for the remote URL.
    """
    try:
        tbx = pysam.TabixFile(url)
        return tbx
    except OSError as e:
        print(f"\n[Error] Could not connect to EBI: {e}")
        return None


def process_bed_regions(input_file, tissue_name, output_file, padding):
    if tissue_name not in TISSUE_MAP:
        print(f"Error: Unknown tissue '{tissue_name}'.")
        print(f"Available keys: {list(TISSUE_MAP.keys())}")
        sys.exit(1)

    ebi_tissue_id = TISSUE_MAP[tissue_name]
    file_url = get_dataset_url(ebi_tissue_id)

    print(f"Connecting to {file_url}...")
    tbx = get_remote_tabix_handle(file_url)
    if tbx is None:
        sys.exit(1)

    print(f"Loading regions from {input_file}...")
    try:
        regions = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
    except Exception as e:
        print(f"Error reading BED file: {e}")
        sys.exit(1)

    final_results = []
    print(f"Querying {len(regions)} regions with +/- {padding}bp padding...")

    count = 0
    total_snps = 0

    for _, row in regions.iterrows():
        chrom = str(row['chrom']).replace('chr', '')
        start = int(row['start'])
        end = int(row['end'])

        q_start = max(1, start - padding)
        q_end = end + padding

        # Retry logic
        max_retries = 3
        fetched_rows = []
        success = False

        for attempt in range(max_retries):
            try:
                # Use CaptureStdErr to catch htslib errors printed to stderr
                with CaptureStdErr() as c:
                    fetched_rows = list(tbx.fetch(chrom, q_start, q_end))

                # 1. Check for stderr errors (Explicit C-level errors)
                if c.output and ("[E::" in c.output or "Illegal seek" in c.output or "Failed" in c.output):
                    raise RuntimeError(f"htslib error: {c.output.strip()}")

                # 2. Check for Empty Result (Potential Silent Failure)
                if not fetched_rows:
                    if attempt < max_retries - 1:
                        raise RuntimeError("Empty result (potential silent network drop). Retrying...")
                    else:
                        # If still empty on last attempt, accept it as truly empty
                        success = True
                        break

                success = True
                break

            except ValueError:
                # Region likely not in file (e.g., chrM), don't retry
                fetched_rows = []
                success = True
                break
            except Exception as e:
                # Print short warning
                print(f"\n  [Retry {attempt + 1}/{max_retries}] {chrom}:{q_start}-{q_end}: {e}")

                # Cleanup and Reconnect
                try:
                    tbx.close()
                except:
                    pass

                time.sleep(2)  # Wait 2 seconds before reconnecting
                tbx = get_remote_tabix_handle(file_url)
                if tbx is None:
                    print("  > Reconnection failed.")
                    break

        if not success:
            print(f"  > Skipping region {chrom}:{q_start}-{q_end} after failures.")
            continue

        input_region_str = f"chr{chrom}:{start}-{end}"

        for line in fetched_rows:
            cols = line.split('\t')
            try:
                se_val = cols[COL_IDX['se']]
                varbeta = float(se_val) ** 2 if se_val != 'NA' else None
                maf_val = cols[COL_IDX['maf']]
                maf = float(maf_val) if maf_val != 'NA' else None
                an_val = cols[COL_IDX['an']]
                N = float(an_val) / 2 if (an_val != 'NA') else 0

                res = {
                    'input_region': input_region_str,
                    'gtex_snp_id': f"chr{cols[COL_IDX['chromosome']]}:{cols[COL_IDX['position']]}",
                    'target_gene_id': cols[COL_IDX['molecular_trait_id']],
                    'beta': float(cols[COL_IDX['beta']]),
                    'varbeta': varbeta,
                    'pval': float(cols[COL_IDX['pvalue']]),
                    'maf': maf,
                    'N': N
                }
                final_results.append(res)
                total_snps += 1
            except (ValueError, IndexError):
                continue

        count += 1
        if count % 10 == 0:
            print(f"Processed {count}/{len(regions)} regions | Found {total_snps} SNPs...", end='\r')

    if final_results:
        print(f"\nProcessing complete. Saving {total_snps} SNPs to {output_file}...")
        df = pd.DataFrame(final_results)
        df.to_csv(output_file, index=False)
        print(f"Success!")
    else:
        print("\nNo data found matching the input regions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fetch GTEx Stats via Remote Tabix')
    parser.add_argument('input_file', help='Input BED file (chrom, start, end)')
    parser.add_argument('--tissue', '-t', required=True, help='Tissue Key (e.g., Adipocytes)')
    parser.add_argument('--output', '-o', default='gtex_tabix_data.csv')
    parser.add_argument('--padding', '-p', type=int, default=1000)

    args = parser.parse_args()
    process_bed_regions(args.input_file, args.tissue, args.output, args.padding)