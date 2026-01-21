library(Biostrings)
library(ggseqlogo)
library(ggplot2)

# Read FASTA file using Biostrings
#NHDF
m6A_all_NHDF_filt <- readRNAStringSet("a.NHDF_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20read.10stoich.expanded.getfasta.allupper.U.fasta")
m6A_DRACH_NHDF_filt <- readRNAStringSet("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif.20read.10stoich.expanded.getfasta.allupper.U.fasta")
#HD10.6
m6A_all_HD10_6_filt <- readRNAStringSet("a.HD10_6_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20read.10stoich.expanded.getfasta.allupper.U.fasta")
m6A_DRACH_HD10_6_filt <- readRNAStringSet("HD10_6_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif.20read.10stoich.expanded.getfasta.allupper.U.fasta")


# Convert to character vector
#NHDF
m6A_all_NHDF_filt_seq <- as.character(m6A_all_NHDF_filt)
m6A_DRACH_NHDF_filt_seq <- as.character(m6A_DRACH_NHDF_filt)
#HD10.6
m6A_all_HD10_6_filt_seq <- as.character(m6A_all_HD10_6_filt)
m6A_DRACH_HD10_6_filt_seq <- as.character(m6A_DRACH_HD10_6_filt)



# Plots
#NHDF
logo_all_NHDF_filt <- ggseqlogo(m6A_all_NHDF_filt_seq, method="prob") +
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
logo_all_NHDF_filt

logo_DRACH_NHDF_filt <- ggseqlogo(m6A_DRACH_NHDF_filt_seq, method="prob") +
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
logo_DRACH_NHDF_filt


#HD10.6
logo_all_HD_filt <- ggseqlogo(m6A_all_HD10_6_filt_seq, method="prob") +
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
logo_all_HD_filt

logo_DRACH_HD_filt <- ggseqlogo(m6A_DRACH_HD10_6_filt_seq, method="prob") +
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
logo_DRACH_HD_filt



#Saving plots

pdf(file = "D) m6A_all in NHDF_DMSO in hg38,0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_all_NHDF_filt
dev.off()

pdf(file = "D) m6A_DRACH in NHDF_DMSO in hg38,0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_DRACH_NHDF_filt
dev.off()


pdf(file = "F) m6A_all in HD10_6_DMSO in hg38,0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_all_HD_filt
dev.off()

pdf(file = "F) m6A_DRACH in HD10_6_DMSO in hg38,0.9.0,filtered,20reads10stoich,motif.pdf", width = 10, height = 5)
logo_DRACH_HD_filt
dev.off()


