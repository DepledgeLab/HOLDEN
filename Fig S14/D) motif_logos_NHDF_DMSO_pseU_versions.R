library(Biostrings)
library(ggseqlogo)
library(ggplot2)

# Read FASTA file using Biostrings
#Versions
pseU_NHDF_DMSO_hg38_7_filt <- readRNAStringSet("unique_7_filtered.expanded.getfasta.allupper.U.fasta")
pseU_NHDF_DMSO_hg38_8_9_filt <- readRNAStringSet("unique_8_9_filtered.expanded.getfasta.allupper.U.fasta")
pseU_NHDF_DMSO_hg38_7_8_9_filt <- readRNAStringSet("7_8_9_filtered.expanded.getfasta.allupper.U.fasta")


# Convert to character vector
pseU_NHDF_DMSO_hg38_7_filt_seq <- as.character(pseU_NHDF_DMSO_hg38_7_filt)
pseU_NHDF_DMSO_hg38_8_9_filt_seq <- as.character(pseU_NHDF_DMSO_hg38_8_9_filt)
pseU_NHDF_DMSO_hg38_7_8_9_filt_seq <- as.character(pseU_NHDF_DMSO_hg38_7_8_9_filt)

# Plots

logo7 <- ggseqlogo(pseU_NHDF_DMSO_hg38_7_filt_seq, method="prob") +
  scale_x_continuous(
    breaks = 1:5,           # Original positions (1 to 5)
    labels = c(-2, -1, 0, 1, 2)  # Your desired labels
  ) +
  xlab("Position") +
  ylab("Frequency") +
  theme(
    axis.title = element_text(size = 26),   # Axis titles
    axis.text = element_text(size = 26)     # Axis tick labels
  )
logo7


logo8_9 <- ggseqlogo(pseU_NHDF_DMSO_hg38_8_9_filt_seq, method="prob") +
  scale_x_continuous(
    breaks = 1:5,           # Original positions (1 to 5)
    labels = c(-2, -1, 0, 1, 2)  # Your desired labels
  ) +
  xlab("Position")+
  ylab("Frequency") +
  theme(
    axis.title = element_text(size = 26),   # Axis titles
    axis.text = element_text(size = 26)     # Axis tick labels
  )
logo8_9


logo7_8_9 <- ggseqlogo(pseU_NHDF_DMSO_hg38_7_8_9_filt_seq, method="prob") +
  scale_x_continuous(
    breaks = 1:5,           # Original positions (1 to 5)
    labels = c(-2, -1, 0, 1, 2)  # Your desired labels
  ) +
  xlab("Position")+
  ylab("Frequency") +
  theme(
    axis.title = element_text(size = 26),   # Axis titles
    axis.text = element_text(size = 26)     # Axis tick labels
  )
logo7_8_9


#Saving plots

pdf(file = "D) pseU in NHDF_DMSO in hg38,0.7.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo7
dev.off()

pdf(file = "D) pseU in NHDF_DMSO in hg38,0.8.0_0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo8_9
dev.off()

pdf(file = "D) pseU in NHDF_DMSO in hg38,0.7.0_0.8.0_0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo7_8_9
dev.off()
