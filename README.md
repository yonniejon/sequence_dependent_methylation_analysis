# Sequence Dependent Methylation Analysis

This repository contains tools and scripts for analyzing sequence-dependent DNA methylation patterns using whole-genome bisulfite sequencing (WGBS) data.

## Overview

This pipeline processes aligned methylation data (BAM) to extract sequence features and correlate them with observed fragment-level methylation patterns.

---

## Data Sources

The analysis is based on the human DNA methylation atlas.

* **Source:** [Loyfer et al., Nature 2023](https://www.nature.com/articles/s41586-022-05580-6)
* **Dataset:** The BAM files used in this analysis are obtained from the **Loyfer Atlas** via the European Genome-phenome Archive (EGA).
* **EGA Study Accession:** [EGAS00001006791](https://ega-archive.org/studies/EGAS00001006791)

---

## Usage Instructions

All processing steps and command-line instructions are documented in the file `execute_pipeline.txt`. To reproduce the analysis or run it on new data, please follow the sequence of commands provided in that file.

### Steps for Execution:

1. **Clone the Repository**
   Use `git clone https://github.com/yonniejon/sequence_dependent_methylation_analysis.git` to download the project.

2. **Data Preparation**
   * Download the required BAM files from the EGA study EGAS00001006791.
   * Ensure they are placed in the expected directory structure as defined in the scripts.

3. **Run the Pipeline**
   * Open the file **execute_pipeline.txt**.
   * Execute the steps sequentially. This typically includes:
     * Pre-processing and filtering of BAM files.
     * Extracting methylation states at specific CpG sites.
     * Calculating sequence-dependent features.
     * Downstream statistical analysis.

---

## Citation

If you use this code or the Loyfer atlas data, please cite:

Rosenski, J., Sabag, O., et al. **The genetic basis for DNA methylation variation across tissues and development.** (2025).  https://doi.org/10.1101/2025.09.15.675351
