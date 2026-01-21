library(data.table)
library(dplyr)
library(tidyr)
library (ggplot2)
library(scales)
library(patchwork)


COV <- 20
FREQ0 <- 0
FREQ5 <- 5
FREQ10 <- 10
FREQ50 <- 50



### Read in transcriptome summary
base <- fread("transcript_regions_summary_gencode_v47.tsv", header = TRUE)
base$result <- NA


### Function to process, all stoichiometries
process_modification_data0 <- function(modification_file, base_data, COV, FREQ0) {
  # Read in modification data
  mod_data <- fread(modification_file, header = FALSE)
  colnames(mod_data) <- c("chrom", "start", "end", "mod", "score", "strand", 
                          "start_tmp", "end_tmp", "color", "valid_cov", "frac_mod", 
                          "N_mod", "N_canon", "N_other", "N_del", "N_fail", "N_diff", "N_nocall")
  
  # Clean chromosome names
  mod_data$chrom <- sub("\\|.*", "", mod_data$chrom)
  
  # Filter data based on coverage and fraction of modification
  mod_data_filtered <- mod_data %>%
    filter(valid_cov >= COV, frac_mod > FREQ0)
  
  # Perform the join and calculate results
  mod_data_processed <- mod_data_filtered %>%
    left_join(base_data, by = c("chrom" = "TranscriptID")) %>%
    mutate(result = case_when(
      is.na(start) ~ NA_real_,
      start < utr5_size ~ start / utr5_size,
      start < (utr5_size + cds_size) ~ (start - utr5_size) / cds_size + 1,
      TRUE ~ (start - utr5_size - cds_size) / utr3_size + 2
    ))
  
  # Filter again after calculation
  mod_data_processed$frac_mod <- as.numeric(mod_data_processed$frac_mod)
  mod_data_processed$valid_cov <- as.numeric(mod_data_processed$valid_cov)
  mod_data_parsed <- mod_data_processed %>%
    filter(valid_cov >= COV & frac_mod > FREQ0)
  
  return(mod_data_parsed)
}

### Function to process, less or equal to 5 stoichiometry
process_modification_data5 <- function(modification_file, base_data, COV, FREQ5) {
  # Read in modification data
  mod_data <- fread(modification_file, header = FALSE)
  colnames(mod_data) <- c("chrom", "start", "end", "mod", "score", "strand", 
                          "start_tmp", "end_tmp", "color", "valid_cov", "frac_mod", 
                          "N_mod", "N_canon", "N_other", "N_del", "N_fail", "N_diff", "N_nocall")
  
  # Clean chromosome names
  mod_data$chrom <- sub("\\|.*", "", mod_data$chrom)
  
  # Filter data based on coverage and fraction of modification
  mod_data_filtered <- mod_data %>%
    filter(valid_cov >= COV, frac_mod <= FREQ5)
  
  # Perform the join and calculate results
  mod_data_processed <- mod_data_filtered %>%
    left_join(base_data, by = c("chrom" = "TranscriptID")) %>%
    mutate(result = case_when(
      is.na(start) ~ NA_real_,
      start < utr5_size ~ start / utr5_size,
      start < (utr5_size + cds_size) ~ (start - utr5_size) / cds_size + 1,
      TRUE ~ (start - utr5_size - cds_size) / utr3_size + 2
    ))
  
  # Filter again after calculation
  mod_data_processed$frac_mod <- as.numeric(mod_data_processed$frac_mod)
  mod_data_processed$valid_cov <- as.numeric(mod_data_processed$valid_cov)
  mod_data_parsed <- mod_data_processed %>%
    filter(valid_cov >= COV & frac_mod <= FREQ5)
  
  return(mod_data_parsed)
}

### Function to process, less than 10 stoichiometry
process_modification_data_less_10 <- function(modification_file, base_data, COV, FREQ10) {
  # Read in modification data
  mod_data <- fread(modification_file, header = FALSE)
  colnames(mod_data) <- c("chrom", "start", "end", "mod", "score", "strand", 
                          "start_tmp", "end_tmp", "color", "valid_cov", "frac_mod", 
                          "N_mod", "N_canon", "N_other", "N_del", "N_fail", "N_diff", "N_nocall")
  
  # Clean chromosome names
  mod_data$chrom <- sub("\\|.*", "", mod_data$chrom)
  
  # Filter data based on coverage and fraction of modification
  mod_data_filtered <- mod_data %>%
    filter(valid_cov >= COV, frac_mod < FREQ10)
  
  # Perform the join and calculate results
  mod_data_processed <- mod_data_filtered %>%
    left_join(base_data, by = c("chrom" = "TranscriptID")) %>%
    mutate(result = case_when(
      is.na(start) ~ NA_real_,
      start < utr5_size ~ start / utr5_size,
      start < (utr5_size + cds_size) ~ (start - utr5_size) / cds_size + 1,
      TRUE ~ (start - utr5_size - cds_size) / utr3_size + 2
    ))
  
  # Filter again after calculation
  mod_data_processed$frac_mod <- as.numeric(mod_data_processed$frac_mod)
  mod_data_processed$valid_cov <- as.numeric(mod_data_processed$valid_cov)
  mod_data_parsed <- mod_data_processed %>%
    filter(valid_cov >= COV & frac_mod < FREQ10)
  
  return(mod_data_parsed)
}

### Function to process, more than or equal to 10 stoichiometry
process_modification_data_more_10 <- function(modification_file, base_data, COV, FREQ10) {
  # Read in modification data
  mod_data <- fread(modification_file, header = FALSE)
  colnames(mod_data) <- c("chrom", "start", "end", "mod", "score", "strand", 
                          "start_tmp", "end_tmp", "color", "valid_cov", "frac_mod", 
                          "N_mod", "N_canon", "N_other", "N_del", "N_fail", "N_diff", "N_nocall")
  
  # Clean chromosome names
  mod_data$chrom <- sub("\\|.*", "", mod_data$chrom)
  
  # Filter data based on coverage and fraction of modification
  mod_data_filtered <- mod_data %>%
    filter(valid_cov >= COV, frac_mod >= FREQ10)
  
  # Perform the join and calculate results
  mod_data_processed <- mod_data_filtered %>%
    left_join(base_data, by = c("chrom" = "TranscriptID")) %>%
    mutate(result = case_when(
      is.na(start) ~ NA_real_,
      start < utr5_size ~ start / utr5_size,
      start < (utr5_size + cds_size) ~ (start - utr5_size) / cds_size + 1,
      TRUE ~ (start - utr5_size - cds_size) / utr3_size + 2
    ))
  
  # Filter again after calculation
  mod_data_processed$frac_mod <- as.numeric(mod_data_processed$frac_mod)
  mod_data_processed$valid_cov <- as.numeric(mod_data_processed$valid_cov)
  mod_data_parsed <- mod_data_processed %>%
    filter(valid_cov >= COV & frac_mod >= FREQ10)
  
  return(mod_data_parsed)
}

### Function to process, more than or equal to 50 stoichiometry
process_modification_data_50 <- function(modification_file, base_data, COV, FREQ50) {
  # Read in modification data
  mod_data <- fread(modification_file, header = FALSE)
  colnames(mod_data) <- c("chrom", "start", "end", "mod", "score", "strand", 
                          "start_tmp", "end_tmp", "color", "valid_cov", "frac_mod", 
                          "N_mod", "N_canon", "N_other", "N_del", "N_fail", "N_diff", "N_nocall")
  
  # Clean chromosome names
  mod_data$chrom <- sub("\\|.*", "", mod_data$chrom)
  
  # Filter data based on coverage and fraction of modification
  mod_data_filtered <- mod_data %>%
    filter(valid_cov >= COV, frac_mod >= FREQ50)
  
  # Perform the join and calculate results
  mod_data_processed <- mod_data_filtered %>%
    left_join(base_data, by = c("chrom" = "TranscriptID")) %>%
    mutate(result = case_when(
      is.na(start) ~ NA_real_,
      start < utr5_size ~ start / utr5_size,
      start < (utr5_size + cds_size) ~ (start - utr5_size) / cds_size + 1,
      TRUE ~ (start - utr5_size - cds_size) / utr3_size + 2
    ))
  
  # Filter again after calculation
  mod_data_processed$frac_mod <- as.numeric(mod_data_processed$frac_mod)
  mod_data_processed$valid_cov <- as.numeric(mod_data_processed$valid_cov)
  mod_data_parsed <- mod_data_processed %>%
    filter(valid_cov >= COV & frac_mod >= FREQ50)
  
  return(mod_data_parsed)
}



#####generating files with different stoichiometries
NHDF_filtered_all_DMSO_parsed_0 <- process_modification_data0("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ0)
NHDF_filtered_all_DMSO_parsed_5 <- process_modification_data5("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ5)
NHDF_filtered_all_DMSO_parsed_less10 <- process_modification_data_less_10("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ10)
NHDF_filtered_all_DMSO_parsed_more10 <- process_modification_data_more_10("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ10)
NHDF_filtered_all_DMSO_parsed_50 <- process_modification_data_50("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ50)

HD_filtered_all_DMSO_parsed_0 <- process_modification_data0("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ0)
HD_filtered_all_DMSO_parsed_5 <- process_modification_data5("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ5)
HD_filtered_all_DMSO_parsed_less10 <- process_modification_data_less_10("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ10)
HD_filtered_all_DMSO_parsed_more10 <- process_modification_data_more_10("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ10)
HD_filtered_all_DMSO_parsed_50 <- process_modification_data_50("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.noFilt.motif.stoich0.bed", base, COV, FREQ50)



# Rescale UTR and CDS
#NHDF
utr5.SF_N_DMSO_0 <- median(NHDF_filtered_all_DMSO_parsed_0$utr5_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_0$cds_size, na.rm = TRUE)
utr3.SF_N_DMSO_0 <- median(NHDF_filtered_all_DMSO_parsed_0$utr3_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_0$cds_size, na.rm = TRUE)

utr5.SF_N_DMSO_5 <- median(NHDF_filtered_all_DMSO_parsed_5$utr5_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_5$cds_size, na.rm = TRUE)
utr3.SF_N_DMSO_5 <- median(NHDF_filtered_all_DMSO_parsed_5$utr3_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_5$cds_size, na.rm = TRUE)

utr5.SF_N_DMSO_l10 <- median(NHDF_filtered_all_DMSO_parsed_less10$utr5_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_less10$cds_size, na.rm = TRUE)
utr3.SF_N_DMSO_l10 <- median(NHDF_filtered_all_DMSO_parsed_less10$utr3_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_less10$cds_size, na.rm = TRUE)

utr5.SF_N_DMSO_m10 <- median(NHDF_filtered_all_DMSO_parsed_more10$utr5_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_more10$cds_size, na.rm = TRUE)
utr3.SF_N_DMSO_m10 <- median(NHDF_filtered_all_DMSO_parsed_more10$utr3_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_more10$cds_size, na.rm = TRUE)

utr5.SF_N_DMSO_50 <- median(NHDF_filtered_all_DMSO_parsed_50$utr5_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_50$cds_size, na.rm = TRUE)
utr3.SF_N_DMSO_50 <- median(NHDF_filtered_all_DMSO_parsed_50$utr3_size, na.rm = TRUE) / median(NHDF_filtered_all_DMSO_parsed_50$cds_size, na.rm = TRUE)

#HD10.6
utr5.SF_H_DMSO_0 <- median(HD_filtered_all_DMSO_parsed_0$utr5_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_0$cds_size, na.rm = TRUE)
utr3.SF_H_DMSO_0 <- median(HD_filtered_all_DMSO_parsed_0$utr3_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_0$cds_size, na.rm = TRUE)

utr5.SF_H_DMSO_5 <- median(HD_filtered_all_DMSO_parsed_5$utr5_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_5$cds_size, na.rm = TRUE)
utr3.SF_H_DMSO_5 <- median(HD_filtered_all_DMSO_parsed_5$utr3_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_5$cds_size, na.rm = TRUE)

utr5.SF_H_DMSO_l10 <- median(HD_filtered_all_DMSO_parsed_less10$utr5_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_less10$cds_size, na.rm = TRUE)
utr3.SF_H_DMSO_l10 <- median(HD_filtered_all_DMSO_parsed_less10$utr3_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_less10$cds_size, na.rm = TRUE)

utr5.SF_H_DMSO_m10 <- median(HD_filtered_all_DMSO_parsed_more10$utr5_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_more10$cds_size, na.rm = TRUE)
utr3.SF_H_DMSO_m10 <- median(HD_filtered_all_DMSO_parsed_more10$utr3_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_more10$cds_size, na.rm = TRUE)

utr5.SF_H_DMSO_50 <- median(HD_filtered_all_DMSO_parsed_50$utr5_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_50$cds_size, na.rm = TRUE)
utr3.SF_H_DMSO_50 <- median(HD_filtered_all_DMSO_parsed_50$utr3_size, na.rm = TRUE) / median(HD_filtered_all_DMSO_parsed_50$cds_size, na.rm = TRUE)


# Split data based on result
#NHDF
utr5.N_DMSO.dist_0 <- NHDF_filtered_all_DMSO_parsed_0[NHDF_filtered_all_DMSO_parsed_0$result < 1, ]
cds.N_DMSO.dist_0 <- NHDF_filtered_all_DMSO_parsed_0[NHDF_filtered_all_DMSO_parsed_0$result < 2 & NHDF_filtered_all_DMSO_parsed_0$result >= 1, ]
utr3.N_DMSO.dist_0 <- NHDF_filtered_all_DMSO_parsed_0[NHDF_filtered_all_DMSO_parsed_0$result >= 2, ]

utr5.N_DMSO.dist_5 <- NHDF_filtered_all_DMSO_parsed_5[NHDF_filtered_all_DMSO_parsed_5$result < 1, ]
cds.N_DMSO.dist_5 <- NHDF_filtered_all_DMSO_parsed_5[NHDF_filtered_all_DMSO_parsed_5$result < 2 & NHDF_filtered_all_DMSO_parsed_5$result >= 1, ]
utr3.N_DMSO.dist_5 <- NHDF_filtered_all_DMSO_parsed_5[NHDF_filtered_all_DMSO_parsed_5$result >= 2, ]

utr5.N_DMSO.dist_l10 <- NHDF_filtered_all_DMSO_parsed_less10[NHDF_filtered_all_DMSO_parsed_less10$result < 1, ]
cds.N_DMSO.dist_l10 <- NHDF_filtered_all_DMSO_parsed_less10[NHDF_filtered_all_DMSO_parsed_less10$result < 2 & NHDF_filtered_all_DMSO_parsed_less10$result >= 1, ]
utr3.N_DMSO.dist_l10 <- NHDF_filtered_all_DMSO_parsed_less10[NHDF_filtered_all_DMSO_parsed_less10$result >= 2, ]

utr5.N_DMSO.dist_m10 <- NHDF_filtered_all_DMSO_parsed_more10[NHDF_filtered_all_DMSO_parsed_more10$result < 1, ]
cds.N_DMSO.dist_m10 <- NHDF_filtered_all_DMSO_parsed_more10[NHDF_filtered_all_DMSO_parsed_more10$result < 2 & NHDF_filtered_all_DMSO_parsed_more10$result >= 1, ]
utr3.N_DMSO.dist_m10 <- NHDF_filtered_all_DMSO_parsed_more10[NHDF_filtered_all_DMSO_parsed_more10$result >= 2, ]

utr5.N_DMSO.dist_50 <- NHDF_filtered_all_DMSO_parsed_50[NHDF_filtered_all_DMSO_parsed_50$result < 1, ]
cds.N_DMSO.dist_50 <- NHDF_filtered_all_DMSO_parsed_50[NHDF_filtered_all_DMSO_parsed_50$result < 2 & NHDF_filtered_all_DMSO_parsed_50$result >= 1, ]
utr3.N_DMSO.dist_50 <- NHDF_filtered_all_DMSO_parsed_50[NHDF_filtered_all_DMSO_parsed_50$result >= 2, ]


#HD10.6
utr5.H_DMSO.dist_0 <- HD_filtered_all_DMSO_parsed_0[HD_filtered_all_DMSO_parsed_0$result < 1, ]
cds.H_DMSO.dist_0 <- HD_filtered_all_DMSO_parsed_0[HD_filtered_all_DMSO_parsed_0$result < 2 & HD_filtered_all_DMSO_parsed_0$result >= 1, ]
utr3.H_DMSO.dist_0 <- HD_filtered_all_DMSO_parsed_0[HD_filtered_all_DMSO_parsed_0$result >= 2, ]

utr5.H_DMSO.dist_5 <- HD_filtered_all_DMSO_parsed_5[HD_filtered_all_DMSO_parsed_5$result < 1, ]
cds.H_DMSO.dist_5 <- HD_filtered_all_DMSO_parsed_5[HD_filtered_all_DMSO_parsed_5$result < 2 & HD_filtered_all_DMSO_parsed_5$result >= 1, ]
utr3.H_DMSO.dist_5 <- HD_filtered_all_DMSO_parsed_5[HD_filtered_all_DMSO_parsed_5$result >= 2, ]

utr5.H_DMSO.dist_l10 <- HD_filtered_all_DMSO_parsed_less10[HD_filtered_all_DMSO_parsed_less10$result < 1, ]
cds.H_DMSO.dist_l10 <- HD_filtered_all_DMSO_parsed_less10[HD_filtered_all_DMSO_parsed_less10$result < 2 & HD_filtered_all_DMSO_parsed_less10$result >= 1, ]
utr3.H_DMSO.dist_l10 <- HD_filtered_all_DMSO_parsed_less10[HD_filtered_all_DMSO_parsed_less10$result >= 2, ]

utr5.H_DMSO.dist_m10 <- HD_filtered_all_DMSO_parsed_more10[HD_filtered_all_DMSO_parsed_more10$result < 1, ]
cds.H_DMSO.dist_m10 <- HD_filtered_all_DMSO_parsed_more10[HD_filtered_all_DMSO_parsed_more10$result < 2 & HD_filtered_all_DMSO_parsed_more10$result >= 1, ]
utr3.H_DMSO.dist_m10 <- HD_filtered_all_DMSO_parsed_more10[HD_filtered_all_DMSO_parsed_more10$result >= 2, ]

utr5.H_DMSO.dist_50 <- HD_filtered_all_DMSO_parsed_50[HD_filtered_all_DMSO_parsed_50$result < 1, ]
cds.H_DMSO.dist_50 <- HD_filtered_all_DMSO_parsed_50[HD_filtered_all_DMSO_parsed_50$result < 2 & HD_filtered_all_DMSO_parsed_50$result >= 1, ]
utr3.H_DMSO.dist_50 <- HD_filtered_all_DMSO_parsed_50[HD_filtered_all_DMSO_parsed_50$result >= 2, ]



# Rescale
#NHDF
utr5.N_DMSO.dist_0$result <- rescale(utr5.N_DMSO.dist_0$result, to = c(1 - utr5.SF_N_DMSO_0, 1), from = c(0, 1))
utr3.N_DMSO.dist_0$result <- rescale(utr3.N_DMSO.dist_0$result, to = c(2, 2 + utr3.SF_N_DMSO_0), from = c(2, 3))

utr5.N_DMSO.dist_5$result <- rescale(utr5.N_DMSO.dist_5$result, to = c(1 - utr5.SF_N_DMSO_5, 1), from = c(0, 1))
utr3.N_DMSO.dist_5$result <- rescale(utr3.N_DMSO.dist_5$result, to = c(2, 2 + utr3.SF_N_DMSO_5), from = c(2, 3))

utr5.N_DMSO.dist_l10$result <- rescale(utr5.N_DMSO.dist_l10$result, to = c(1 - utr5.SF_N_DMSO_l10, 1), from = c(0, 1))
utr3.N_DMSO.dist_l10$result <- rescale(utr3.N_DMSO.dist_l10$result, to = c(2, 2 + utr3.SF_N_DMSO_l10), from = c(2, 3))

utr5.N_DMSO.dist_m10$result <- rescale(utr5.N_DMSO.dist_m10$result, to = c(1 - utr5.SF_N_DMSO_m10, 1), from = c(0, 1))
utr3.N_DMSO.dist_m10$result <- rescale(utr3.N_DMSO.dist_m10$result, to = c(2, 2 + utr3.SF_N_DMSO_m10), from = c(2, 3))

utr5.N_DMSO.dist_50$result <- rescale(utr5.N_DMSO.dist_50$result, to = c(1 - utr5.SF_N_DMSO_50, 1), from = c(0, 1))
utr3.N_DMSO.dist_50$result <- rescale(utr3.N_DMSO.dist_50$result, to = c(2, 2 + utr3.SF_N_DMSO_50), from = c(2, 3))



#HD10.6
utr5.H_DMSO.dist_0$result <- rescale(utr5.H_DMSO.dist_0$result, to = c(1 - utr5.SF_H_DMSO_0, 1), from = c(0, 1))
utr3.H_DMSO.dist_0$result <- rescale(utr3.H_DMSO.dist_0$result, to = c(2, 2 + utr3.SF_H_DMSO_0), from = c(2, 3))

utr5.H_DMSO.dist_5$result <- rescale(utr5.H_DMSO.dist_5$result, to = c(1 - utr5.SF_H_DMSO_5, 1), from = c(0, 1))
utr3.H_DMSO.dist_5$result <- rescale(utr3.H_DMSO.dist_5$result, to = c(2, 2 + utr3.SF_H_DMSO_5), from = c(2, 3))

utr5.H_DMSO.dist_l10$result <- rescale(utr5.H_DMSO.dist_l10$result, to = c(1 - utr5.SF_H_DMSO_l10, 1), from = c(0, 1))
utr3.H_DMSO.dist_l10$result <- rescale(utr3.H_DMSO.dist_l10$result, to = c(2, 2 + utr3.SF_H_DMSO_l10), from = c(2, 3))

utr5.H_DMSO.dist_m10$result <- rescale(utr5.H_DMSO.dist_m10$result, to = c(1 - utr5.SF_H_DMSO_m10, 1), from = c(0, 1))
utr3.H_DMSO.dist_m10$result <- rescale(utr3.H_DMSO.dist_m10$result, to = c(2, 2 + utr3.SF_H_DMSO_m10), from = c(2, 3))

utr5.H_DMSO.dist_50$result <- rescale(utr5.H_DMSO.dist_50$result, to = c(1 - utr5.SF_H_DMSO_50, 1), from = c(0, 1))
utr3.H_DMSO.dist_50$result <- rescale(utr3.H_DMSO.dist_50$result, to = c(2, 2 + utr3.SF_H_DMSO_50), from = c(2, 3))



# Combine results
N_DMSO.metagene.coord_0 <- c(utr5.N_DMSO.dist_0$result, cds.N_DMSO.dist_0$result, utr3.N_DMSO.dist_0$result)
N_DMSO.metagene.coord_5 <- c(utr5.N_DMSO.dist_5$result, cds.N_DMSO.dist_5$result, utr3.N_DMSO.dist_5$result)
N_DMSO.metagene.coord_l10 <- c(utr5.N_DMSO.dist_l10$result, cds.N_DMSO.dist_l10$result, utr3.N_DMSO.dist_l10$result)
N_DMSO.metagene.coord_m10 <- c(utr5.N_DMSO.dist_m10$result, cds.N_DMSO.dist_m10$result, utr3.N_DMSO.dist_m10$result)
N_DMSO.metagene.coord_50 <- c(utr5.N_DMSO.dist_50$result, cds.N_DMSO.dist_50$result, utr3.N_DMSO.dist_50$result)


H_DMSO.metagene.coord_0 <- c(utr5.H_DMSO.dist_0$result, cds.H_DMSO.dist_0$result, utr3.H_DMSO.dist_0$result)
H_DMSO.metagene.coord_5 <- c(utr5.H_DMSO.dist_5$result, cds.H_DMSO.dist_5$result, utr3.H_DMSO.dist_5$result)
H_DMSO.metagene.coord_l10 <- c(utr5.H_DMSO.dist_l10$result, cds.H_DMSO.dist_l10$result, utr3.H_DMSO.dist_l10$result)
H_DMSO.metagene.coord_m10 <- c(utr5.H_DMSO.dist_m10$result, cds.H_DMSO.dist_m10$result, utr3.H_DMSO.dist_m10$result)
H_DMSO.metagene.coord_50 <- c(utr5.H_DMSO.dist_50$result, cds.H_DMSO.dist_50$result, utr3.H_DMSO.dist_50$result)



#for plots
m6A_all <- c(N_DMSO.metagene.coord_0, H_DMSO.metagene.coord_0)
Cells_all <- c(rep("NHDF", length(N_DMSO.metagene.coord_0)), 
               rep("HD10.6", length(H_DMSO.metagene.coord_0))) 
df_0 <- data.frame(m6A_all, Cells_all)


m6A_less5 <- c(N_DMSO.metagene.coord_5, H_DMSO.metagene.coord_5)
Cells_less5 <- c(rep("NHDF", length(N_DMSO.metagene.coord_5)), 
                 rep("HD10.6", length(H_DMSO.metagene.coord_5))) 
df_5 <- data.frame(m6A_less5, Cells_less5)


m6A_less10 <- c(N_DMSO.metagene.coord_l10, H_DMSO.metagene.coord_l10)
Cells_less10 <- c(rep("NHDF", length(N_DMSO.metagene.coord_l10)), 
                  rep("HD10.6", length(H_DMSO.metagene.coord_l10))) 
df_l10 <- data.frame(m6A_less10, Cells_less10)


m6A_more10 <- c(N_DMSO.metagene.coord_m10, H_DMSO.metagene.coord_m10)
Cells_more10 <- c(rep("NHDF", length(N_DMSO.metagene.coord_m10)), 
                  rep("HD10.6", length(H_DMSO.metagene.coord_m10))) 
df_m10 <- data.frame(m6A_more10, Cells_more10)


m6A_more50 <- c(N_DMSO.metagene.coord_50, H_DMSO.metagene.coord_50)
Cells_more50 <- c(rep("NHDF", length(N_DMSO.metagene.coord_50)), 
                  rep("HD10.6", length(H_DMSO.metagene.coord_50))) 
df_50 <- data.frame(m6A_more50, Cells_more50)



# Plot
all_stoich <- qplot(m6A_all, data = df_0, geom = "density", color = Cells_all, size = I(1.2)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("NHDF" = "darkred", "HD10.6" = "red")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),  # Axis labels
    axis.text = element_text(size = 14),   # Axis tick labels
    legend.text = element_text(size = 14), # Legend labels
    legend.title = element_text(size = 16) # Legend title
  )
all_stoich


less_5 <- qplot(m6A_less5, data = df_5, geom = "density", color = Cells_less5, size = I(1.2)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("NHDF" = "darkred", "HD10.6" = "red")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),  # Axis labels
    axis.text = element_text(size = 14),   # Axis tick labels
    legend.text = element_text(size = 14), # Legend labels
    legend.title = element_text(size = 16) # Legend title
  )
less_5


less_10 <- qplot(m6A_less10, data = df_l10, geom = "density", color = Cells_less10, size = I(1.2)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("NHDF" = "darkred", "HD10.6" = "red")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),  # Axis labels
    axis.text = element_text(size = 14),   # Axis tick labels
    legend.text = element_text(size = 14), # Legend labels
    legend.title = element_text(size = 16) # Legend title
  )
less_10


more_10 <- qplot(m6A_more10, data = df_m10, geom = "density", color = Cells_more10, size = I(1.2)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("NHDF" = "darkred", "HD10.6" = "red")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),  # Axis labels
    axis.text = element_text(size = 14),   # Axis tick labels
    legend.text = element_text(size = 14), # Legend labels
    legend.title = element_text(size = 16) # Legend title
  )
more_10


more_50 <- qplot(m6A_more50, data = df_50, geom = "density", color = Cells_more50, size = I(1.2)) + 
  scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("NHDF" = "darkred", "HD10.6" = "red")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),  # Axis labels
    axis.text = element_text(size = 14),   # Axis tick labels
    legend.text = element_text(size = 14), # Legend labels
    legend.title = element_text(size = 16) # Legend title
  )
more_50


combined_plot <- all_stoich + more_50 + more_10 + less_10 + less_5 + plot_layout(ncol = 1)
combined_plot


pdf(file = "Unfiltered_m6A_DRACH_NHDF_HD10.6_stoichiometry_differences_tx_morethan20reads_all_supplementary.pdf", width = 9, height = 25)
combined_plot
dev.off()
