#!/usr/bin/env Rscript

# Load required libraries
if (!require("coloc")) install.packages("coloc")
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")
if (!require("R.utils")) install.packages("R.utils")

library(coloc)
library(data.table)
library(dplyr)
library(R.utils)

# ==============================================================================
# CONFIGURATION & METADATA
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_coloc.R <sd_asm_file> <gwas_file> <output_file> [N_SD_ASM_FALLBACK]")
}

SD_ASM_FILE <- args[1]
GWAS_FILE   <- args[2]
OUTPUT_FILE <- args[3]
N_SD_ASM_FALLBACK <- if (length(args) >= 4) as.numeric(args[4]) else NA
GNOMAD_REF  <- "/cs/zbio/jrosensk/gnomAD/gnomAD.all.hg38.bed.gz" 
WINDOW_SIZE <- 20000 

study_metadata <- list(
  "22961000-GCST005581" = list(type = "cc", cases = 2861, controls = 8514, trait = "Primary biliary cirrhosis"),
  "23143596-GCST005569" = list(type = "cc", cases = 13838, controls = 33742, trait = "Rheumatoid arthritis (UKB)"),
  "24076602-GCST005531" = list(type = "cc", cases = 14498, controls = 24091, trait = "Multiple sclerosis"),
  "25751624-GCST005536" = list(type = "cc", cases = 6683, controls = 12173, trait = "Type 1 diabetes"),
  "26691988-GCST003219" = list(type = "cc", cases = 16144, controls = 17832, trait = "Age-related macular degeneration"),
  "27416945-GCST003658" = list(type = "quant", n = 16753, trait = "Modified Stumvoll Insulin Sensitivity Index"),
  "28714975-GCST004787" = list(type = "cc", cases = 18467, controls = 45264, trait = "Coronary artery disease"),
  "28887542-GCST005058" = list(type = "quant", n = 9796, trait = "HDL cholesterol"),
  "28887542-GCST005064" = list(type = "quant", n = 9463, trait = "Aspartate aminotransferase levels"),
  "28887542-GCST005066" = list(type = "quant", n = 9803, trait = "Creatinine levels"),
  "28887542-GCST005071" = list(type = "quant", n = 9732, trait = "Insulin-like growth factor 1 levels"),
  "28887542-GCST005072" = list(type = "quant", n = 9818, trait = "Ferritin levels"),
  "29273806-GCST006862" = list(type = "cc", cases = 19954, controls = 107715, trait = "Asthma"),
  "29293525-GCST005353" = list(type = "quant", n = 79, trait = "GLP-1-stimulated insulin secretion"),
  "29358691-GCST005413" = list(type = "cc", cases = 12931, controls = 57196, trait = "Type 2 diabetes"),
  "29404214-GCST006305" = list(type = "quant", n = 2782, trait = "Cholesterol"),
  "29404214-GCST006306" = list(type = "quant", n = 2782, trait = "HDL cholesterol"),
  "32887889-GCST90011813" = list(type = "cc", cases = 762, controls = 410350, trait = "Thyroid cancer"),
  "26831199-GCST003374" = list(type = "cc", cases = 12385, controls = 104780, trait = "Chronic kidney disease"),
  "26831199-GCST003401" = list(type = "quant", n = 118448, trait = "Glomerular filtration rate in non diabetics"),
  "27618448-GCST006227" = list(type = "quant", n = 140886, trait = "Systolic blood pressure"),
  "27618448-GCST006228" = list(type = "quant", n = 140886, trait = "Diastolic blood pressure"),
  "27618448-GCST006230" = list(type = "quant", n = 140886, trait = "Pulse pressure"),
  "27618448-GCST006231" = list(type = "quant", n = 140886, trait = "Mean arterial pressure"),
  "28490609-GCST004487" = list(type = "quant", n = 5567, trait = "Peak insulin response"),
  "28490609-GCST004488" = list(type = "quant", n = 5567, trait = "Insulin secretion rate"),
  "28613276-GCST004734" = list(type = "quant", n = 27850, trait = "Heart rate variability traits"),
  "29212778-GCST005194" = list(type = "cc", cases = 34541, controls = 261984, trait = "Coronary artery disease (UK Biobank)"),
  "GCST90480650" = list(type = "quant", n = 635000, trait = "Platelet count")
)

# ==============================================================================
# 1. GWAS PREPARATION
# ==============================================================================
cat("Loading GWAS data...\n")
gwas <- fread(GWAS_FILE)

file_key <- sub("^([^-]+-[^-]+).*", "\\1", basename(GWAS_FILE))
if (!(file_key %in% names(study_metadata))) {
  file_key <- sub("^([^.]+).*", "\\1", basename(GWAS_FILE))
}
meta <- study_metadata[[file_key]]

# --- 1A. Assign GWAS_TYPE ---
if (!is.null(meta)) {
  GWAS_TYPE <- meta$type
} else {
  has_or_cols <- any(c("hm_odds_ratio", "odds_ratio") %in% names(gwas))
  GWAS_TYPE <- if(has_or_cols) "cc" else "quant"
}

# --- 1B. SAFE COLUMN ASSIGNMENTS ---
# Chromosome
if ("chromosome" %in% names(gwas)) {
  gwas$chrom <- sub("^chr", "", as.character(gwas$chromosome), ignore.case=T)
} else if ("hm_chrom" %in% names(gwas)) {
  gwas$chrom <- sub("^chr", "", as.character(gwas$hm_chrom), ignore.case=T)
} else {
  stop("Chromosome column not found.")
}

# Position
if ("base_pair_location" %in% names(gwas)) {
  gwas$pos <- as.numeric(gwas$base_pair_location)
} else if ("hm_pos" %in% names(gwas)) {
  gwas$pos <- as.numeric(gwas$hm_pos)
} else {
  stop("Position column not found.")
}

# P-value
gwas$p_val <- as.numeric(gwas$p_value)

# EAF/MAF
if ("effect_allele_frequency" %in% names(gwas)) {
  gwas$maf_file <- as.numeric(gwas$effect_allele_frequency)
} else if ("hm_effect_allele_frequency" %in% names(gwas)) {
  gwas$maf_file <- as.numeric(gwas$hm_effect_allele_frequency)
} else {
  gwas$maf_file <- NA_real_
}

# --- 1C. Sample Size ---
n_col_name <- intersect(names(gwas), c("n", "N", "hm_n", "sample_size"))[1]
if (!is.na(n_col_name)) {
  gwas$final_n <- as.numeric(gwas[[n_col_name]])
} else if (!is.null(meta)) {
  gwas$final_n <- if(meta$type == "cc") (meta$cases + meta$controls) else meta$n
} else {
  gwas$final_n <- NA_real_
}

# --- 1D. Case Proportion 's' ---
GWAS_S_VALUE <- NA_real_
if (GWAS_TYPE == "cc") {
  if (all(c("num_cases", "num_controls") %in% names(gwas))) {
    gwas$s_file <- as.numeric(gwas$num_cases) / (as.numeric(gwas$num_cases) + as.numeric(gwas$num_controls))
    GWAS_S_VALUE <- median(gwas$s_file, na.rm=TRUE)
  } else if (all(c("n_cases", "n_controls") %in% names(gwas))) {
    gwas$s_file <- as.numeric(gwas$n_cases) / (as.numeric(gwas$n_cases) + as.numeric(gwas$n_controls))
    GWAS_S_VALUE <- median(gwas$s_file, na.rm=TRUE)
  } else if (!is.null(meta)) {
    GWAS_S_VALUE <- meta$cases / (meta$cases + meta$controls)
  }
}

# ==============================================================================
# 2. MAF, BETA, & VARIANCE CALCULATION
# ==============================================================================
gwas$maf <- gwas$maf_file
if (any(is.na(gwas$maf))) {
  cat("Merging with gnomAD for missing MAF...\n")
  gnomad <- fread(GNOMAD_REF, select = c(1, 3, 7), col.names = c("ref_chrom", "ref_pos", "ref_af"))
  gnomad$ref_chrom <- sub("^chr", "", gnomad$ref_chrom, ignore.case = TRUE)
  gnomad <- gnomad[order(-ref_af), .SD[1], by = .(ref_chrom, ref_pos)]
  gwas <- left_join(gwas, gnomad, by = c("chrom" = "ref_chrom", "pos" = "ref_pos"))
  gwas$maf <- ifelse(is.na(gwas$maf), pmin(gwas$ref_af, 1 - gwas$ref_af), gwas$maf)
  rm(gnomad); gc()
}

# Safe Beta Extraction
if ("beta" %in% names(gwas)) {
  gwas$final_beta <- as.numeric(gwas$beta)
} else if ("hm_beta" %in% names(gwas)) {
  gwas$final_beta <- as.numeric(gwas$hm_beta)
} else if (GWAS_TYPE == "cc" && "odds_ratio" %in% names(gwas)) {
  gwas$final_beta <- log(as.numeric(gwas$odds_ratio))
} else if (GWAS_TYPE == "cc" && "hm_odds_ratio" %in% names(gwas)) {
  gwas$final_beta <- log(as.numeric(gwas$hm_odds_ratio))
} else {
  gwas$final_beta <- NA_real_
}

# Safe SE Extraction
if ("standard_error" %in% names(gwas)) {
  gwas$final_se <- as.numeric(gwas$standard_error)
} else if ("se" %in% names(gwas)) {
  gwas$final_se <- as.numeric(gwas$se)
} else if (!all(is.na(gwas$final_beta)) && "p_value" %in% names(gwas)) {
  gwas$final_se <- abs(gwas$final_beta / qnorm(as.numeric(gwas$p_value)/2, lower.tail=FALSE))
} else {
  gwas$final_se <- NA_real_
}

# Quality filter
gwas <- gwas %>% filter(!is.na(chrom) & !is.na(pos) & !is.na(final_n) & final_n > 1 & !is.na(maf))
if (nrow(gwas) == 0) stop("No usable GWAS SNPs after validation.")

# ==============================================================================
# 3. ANALYSIS (Remaining logic as per source)
# ==============================================================================
cat("Loading SD-ASM data...\n")
sd_asm <- fread(SD_ASM_FILE)
sd_asm$N_sd_final <- if ("N" %in% names(sd_asm)) as.numeric(sd_asm$N) else N_SD_ASM_FALLBACK

gwas$snp_id_match <- paste0("chr", gwas$chrom, ":", gwas$pos)
merged_data <- inner_join(gwas, sd_asm, by = c("snp_id_match" = "snp_id"), suffix = c(".gwas", ".sd")) %>%
  arrange(chrom, pos)

if (nrow(merged_data) == 0) stop("Zero overlapping SNPs found.")

results_list <- list(); i <- 1; total_snps <- nrow(merged_data)
while (i <= total_snps) {
  cur_chrom <- merged_data$chrom[i]; cur_pos <- merged_data$pos[i]
  window_idx <- which(merged_data$chrom == cur_chrom & merged_data$pos >= cur_pos & merged_data$pos <= cur_pos + WINDOW_SIZE)
  
  if (length(window_idx) >= 2) {
    w_data <- merged_data[window_idx, ]
    
    # Newly added to remove dups
    w_data <- w_data %>%
      arrange(meqtl_p_val) %>% 
      distinct(snp_id_match, .keep_all = TRUE) %>%
      arrange(chrom, pos) # Re-sort by position for consistency
    
    # Re-check length after removing duplicates
    if (nrow(w_data) < 2) {
      i <- max(window_idx) + 1
      next
    }
    # end newly added
    
    d1 <- list(snp = w_data$snp_id_match, pvalues = w_data$meqtl_p_val, type = "quant", N = median(w_data$N_sd_final), MAF = w_data$maf, sdY = 1)
    d2 <- list(snp = w_data$snp_id_match, type = GWAS_TYPE, N = median(w_data$final_n), MAF = w_data$maf)
    
    if (all(!is.na(w_data$final_beta)) && all(!is.na(w_data$final_se))) {
      d2$beta <- w_data$final_beta; d2$varbeta <- (w_data$final_se)^2
      if (GWAS_TYPE == "cc") d2$s <- GWAS_S_VALUE else d2$sdY <- 1
    } else {
      d2$pvalues <- w_data$p_val
      if (GWAS_TYPE == "cc") d2$s <- GWAS_S_VALUE else d2$sdY <- 1
    }
    
    res <- tryCatch({ 
      coloc.abf(d1, d2) 
    }, error = function(e) {
      # Print the error message to the console
      cat(paste0("\n[!] Error in region: ", cur_chrom, ":", cur_pos, "\n"))
      cat(paste0("    Message: ", e$message, "\n"))
      return(NULL) 
    })
    
    if (!is.null(res)) {
      top <- res$results[which.max(res$results$SNP.PP.H4), ]
      s <- res$summary
      results_list[[length(results_list) + 1]] <- data.frame(
        region = paste0("chr", cur_chrom, ":", cur_pos, "-", cur_pos + WINDOW_SIZE),
        top_snp = top$snp, top_snp_pp_h4 = top$SNP.PP.H4,
        n_snps = length(window_idx), 
        PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"], PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"]
      )
    }
  }
  i <- max(window_idx) + 1
}

if (length(results_list) > 0) {
  write.csv(bind_rows(results_list) %>% arrange(desc(PP.H4)), OUTPUT_FILE, row.names=FALSE)
  cat("\nProcess complete.\n")
}