#!/usr/bin/env Rscript

# Load required libraries
if (!require("coloc")) install.packages("coloc")
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")

library(coloc)
library(data.table)
library(dplyr)

# ==============================================================================
# CONFIGURATION & CLI PARSING
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("No arguments provided. Using default hardcoded paths.\n")
  cat("Usage: Rscript run_coloc_analysis.R <sd_asm_file> <gtex_file> <output_file> [N_SD_ASM]\n")
  
  SD_ASM_FILE <- "Adipocytes.sd_asm_beta.csv"
  GTEX_FILE   <- "Adipocytes.GTEx.summary_stats_near_sd.csv"
  OUTPUT_FILE <- "Adipocytes_Coloc_Results.csv"
  N_SD_ASM    <- 111 
} else {
  if (length(args) < 3) {
    stop("Error: Please provide at least 3 arguments: <sd_asm_file> <gtex_file> <output_file>")
  }
  SD_ASM_FILE <- args[1]
  GTEX_FILE   <- args[2]
  OUTPUT_FILE <- args[3]
  N_SD_ASM    <- if (length(args) >= 4) as.numeric(args[4]) else 111
}

MIN_OVERLAP_SNPS <- 1  

cat("----------------------------------------------------------\n")
cat("Starting Colocalization Analysis\n")
cat("  SD-ASM File: ", SD_ASM_FILE, "\n")
cat("  GTEx File:   ", GTEX_FILE, "\n")
cat("  Output File: ", OUTPUT_FILE, "\n")
cat("  N (SD-ASM):  ", N_SD_ASM, "\n")
cat("----------------------------------------------------------\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
if (!file.exists(SD_ASM_FILE)) stop(paste("File not found:", SD_ASM_FILE))
if (!file.exists(GTEX_FILE)) stop(paste("File not found:", GTEX_FILE))

cat("Loading SD-ASM data...\n")
sd_asm <- fread(SD_ASM_FILE)
sd_asm$snp_id <- as.character(sd_asm$snp_id)

cat("Loading GTEx data...\n")
gtex <- fread(GTEX_FILE)
gtex$gtex_snp_id <- as.character(gtex$gtex_snp_id)

# ==============================================================================
# PREPARE ANALYSIS LIST
# ==============================================================================
analysis_groups <- gtex %>%
  select(input_region, target_gene_id) %>%
  distinct()

cat("Found", nrow(analysis_groups), "unique Region-Gene pairs to analyze.\n")

results_list <- list()

# ==============================================================================
# RUN COLOC LOOP
# ==============================================================================
for (i in 1:nrow(analysis_groups)) {
  
  region_id <- analysis_groups$input_region[i]
  gene_id   <- analysis_groups$target_gene_id[i]
  
  gtex_subset <- gtex[input_region == region_id & target_gene_id == gene_id]
  merged_data <- inner_join(gtex_subset, sd_asm, by = c("gtex_snp_id" = "snp_id"), suffix = c(".gtex", ".sd"))
  
  if (nrow(merged_data) < MIN_OVERLAP_SNPS) {
    next
  }
  
  # Ensure MAF is valid (coloc crashes if MAF is NA or 0)
  # GTEx usually provides good MAF, but we fill NAs just in case
  merged_data$maf <- as.numeric(merged_data$maf)
  merged_data <- merged_data %>% filter(!is.na(maf) & maf > 0 & maf < 1)
  
  if (nrow(merged_data) < MIN_OVERLAP_SNPS) {
    next # Skip if filtering MAF removed all SNPs
  }

  # 4. Prepare Coloc Datasets
  # Dataset 1: SD-ASM
  # IMPROVEMENT: We use GTEx MAF for SD-ASM dataset too (since SNPs match).
  # This allows coloc to estimate sdY automatically, so we REMOVED 'sdY=1'.
  dataset1 <- list(
    beta = merged_data$beta.sd,
    varbeta = merged_data$varbeta.sd,
    snp = merged_data$gtex_snp_id,
    type = "quant",
    N = N_SD_ASM,
    sdY = 1,
    MAF = merged_data$maf 
  )
  
  # Dataset 2: GTEx
  gtex_N <- median(merged_data$N, na.rm = TRUE)
  if (is.na(gtex_N) || gtex_N < 1) gtex_N <- 500
  
  dataset2 <- list(
    beta = merged_data$beta.gtex,
    varbeta = merged_data$varbeta.gtex,
    snp = merged_data$gtex_snp_id,
    type = "quant",
    N = gtex_N,
    sdY = 1, # GTEx is typically rank-normalized, so sdY=1 is correct here
    MAF = merged_data$maf
  )
  
  res <- tryCatch({
    coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
  }, error = function(e) { 
    cat(sprintf("\n[Coloc Error] Region: %s | Gene: %s | Msg: %s\n", region_id, gene_id, e$message))
    return(NULL) 
  })
  
  if (!is.null(res)) {
    summary_vec <- res$summary
    results_list[[length(results_list) + 1]] <- data.frame(
      input_region = region_id,
      target_gene = gene_id,
      n_snps = nrow(merged_data),
      PP.H0 = summary_vec["PP.H0.abf"],
      PP.H1 = summary_vec["PP.H1.abf"], 
      PP.H2 = summary_vec["PP.H2.abf"], 
      PP.H3 = summary_vec["PP.H3.abf"], 
      PP.H4 = summary_vec["PP.H4.abf"]  
    )
  }
  
  if (i %% 10 == 0) cat(sprintf("Processed %d / %d pairs...\r", i, nrow(analysis_groups)))
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
if (length(results_list) > 0) {
  final_results <- bind_rows(results_list)
  final_results <- final_results %>% arrange(desc(PP.H4))
  write.csv(final_results, OUTPUT_FILE, row.names = FALSE)
  cat("\nSuccess! Results saved to", OUTPUT_FILE, "\n")
  print(head(final_results))
} else {
  cat("\nNo successful colocalization analyses performed.\n")
}
