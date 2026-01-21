library(Biostrings)
library(ggseqlogo)
library(ggplot2)


#NHDF
# Read FASTA file using Biostrings
m6A_all_NHDF_unf <- readRNAStringSet("a.NHDF_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_10stoich.expanded.getfasta.allupper.U.fasta")
m6A_DRACH_NHDF_unf <- readRNAStringSet("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif_20read_10stoich.expanded.getfasta.allupper.U.fasta")
#HD10.6
m6A_all_HD_unf <- readRNAStringSet("a.HD10_6_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_10stoich.expanded.getfasta.allupper.U.fasta")
m6A_DRACH_HD_unf <- readRNAStringSet("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif_20read_10stoich.expanded.getfasta.allupper.U.fasta")


# Convert to character vector
m6A_all_NHDF_unf_seq <- as.character(m6A_all_NHDF_unf)
m6A_DRACH_NHDF_unf_seq <- as.character(m6A_DRACH_NHDF_unf)
#HD10.6
m6A_all_HD_unf_seq <- as.character(m6A_all_HD_unf)
m6A_DRACH_HD_unf_seq <- as.character(m6A_DRACH_HD_unf)


# Plots
#NHDF
logo_all_NHDF_unfilt <- ggseqlogo(m6A_all_NHDF_unf_seq, method="prob") +
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
logo_all_NHDF_unfilt

logo_DRACH_NHDF_unfilt <- ggseqlogo(m6A_DRACH_NHDF_unf_seq, method="prob") +
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
logo_DRACH_NHDF_unfilt


logo_all_HD_unfilt <- ggseqlogo(m6A_all_HD_unf_seq, method="prob") +
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
logo_all_HD_unfilt

logo_DRACH_HD_unfilt <- ggseqlogo(m6A_DRACH_HD_unf_seq, method="prob") +
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
logo_DRACH_HD_unfilt




#Saving plots

pdf(file = "A) m6A_all in NHDF_DMSO in hg38,0.9.0,unfiltered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_all_NHDF_unfilt
dev.off()

pdf(file = "A) m6A_DRACH in NHDF_DMSO in hg38,0.9.0,unfiltered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_DRACH_NHDF_unfilt
dev.off()

pdf(file = "B) m6A_all in HD10_6_DMSO in hg38,0.9.0,unfiltered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_all_HD_unfilt
dev.off()

pdf(file = "B) m6A_DRACH in HD10_6_DMSO in hg38,0.9.0,unfiltered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_DRACH_HD_unfilt
dev.off()


