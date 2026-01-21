library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#datasets
NHDF_DRACH <- na.omit(annot_NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif_sorted_20reads_10stoich_6_dist)
NHDF_a <- na.omit(annot_a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif_sorted_20reads_10stoich_6_dist)
NHDF_DRACH_mono <- na.omit(annot_NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif_sorted_20reads_10stoich_6_dist.monoexonic)
NHDF_a_mono <- na.omit(annot_a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif_sorted_20reads_10stoich_6_dist.monoexonic)

HD_DRACH <- na.omit(annot_HD10_6_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif_sorted_20reads_10stoich_6_dist)
HD_a <- na.omit(annot_a.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif_sorted_20reads_10stoich_6_dist)
HD_DRACH_mono <- na.omit(annot_HD10_6_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif_sorted_20reads_10stoich_6_dist.monoexonic)
HD_a_mono <- na.omit(annot_a.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif_sorted_20reads_10stoich_6_dist.monoexonic)


#dists
NHDF_DRACH_dist <- NHDF_DRACH[NHDF_DRACH$cds_size != 0 & NHDF_DRACH$utr5_size != 0 & NHDF_DRACH$utr3_size != 0, ]
HD_DRACH_dist  <- HD_DRACH[HD_DRACH$cds_size != 0 & HD_DRACH$utr5_size != 0 & HD_DRACH$utr3_size != 0, ]
NHDF_DRACH_mono_dist <- NHDF_DRACH_mono[NHDF_DRACH_mono$cds_size != 0 & NHDF_DRACH_mono$utr5_size != 0 & NHDF_DRACH_mono$utr3_size != 0, ]
HD_DRACH_mono_dist  <- HD_DRACH_mono[HD_DRACH_mono$cds_size != 0 & HD_DRACH_mono$utr5_size != 0 & HD_DRACH_mono$utr3_size != 0, ]

NHDF_a_dist  <- NHDF_a[NHDF_a$cds_size != 0 & NHDF_a$utr5_size != 0 & NHDF_a$utr3_size != 0, ]
HD_a_dist <- HD_a[HD_a$cds_size != 0 & HD_a$utr5_size != 0 & HD_a$utr3_size != 0, ]
NHDF_a_mono_dist  <- NHDF_a_mono[NHDF_a_mono$cds_size != 0 & NHDF_a_mono$utr5_size != 0 & NHDF_a_mono$utr3_size != 0, ]
HD_a_mono_dist  <- HD_a_mono[HD_a_mono$cds_size != 0 & HD_a_mono$utr5_size != 0 & HD_a_mono$utr3_size != 0, ]


# Determine longest length transcript for each gene
#trx_len_N <- NHDF_DRACH_dist$utr5_size + NHDF_DRACH_dist$cds_size + NHDF_DRACH_dist$utr3_size # Determine transcript length
#trx_len_H <- HD_DRACH_dist$utr5_size + HD_DRACH_dist$cds_size + HD_DRACH_dist$utr3_size # Determine transcript length
#trx_len_NHDF_DRACH_mono_dist <- NHDF_DRACH_mono_dist$utr5_size + NHDF_DRACH_mono_dist$cds_size + NHDF_DRACH_mono_dist$utr3_size # Determine transcript length
#trx_len_HD_DRACH_mono_dist <- HD_DRACH_mono_dist$utr5_size + HD_DRACH_mono_dist$cds_size + HD_DRACH_mono_dist$utr3_size # Determine transcript length

#temp_N <- data.frame(NHDF_DRACH_dist$gene_name, NHDF_DRACH_dist$refseqID, trx_len_N)
#temp_H <- data.frame(HD_DRACH_dist$gene_name, HD_DRACH_dist$refseqID, trx_len_H)
#temp_N_mono <- data.frame(NHDF_DRACH_mono_dist$gene_name, NHDF_DRACH_mono_dist$refseqID, trx_len_NHDF_DRACH_mono_dist)
#temp_H_mono <- data.frame(HD_DRACH_mono_dist$gene_name, HD_DRACH_mono_dist$refseqID, trx_len_HD_DRACH_mono_dist)

#colnames(temp_N) <- c("gene_name", "gid", "trx_len_N") 
#colnames(temp_H) <- c("gene_name", "gid", "trx_len_H") 
#colnames(temp_N_mono) <- c("gene_name", "gid", "trx_len_NHDF_DRACH_mono_dist") 
#colnames(temp_H_mono) <- c("gene_name", "gid", "trx_len_HD_DRACH_mono_dist") 

#temp_N <- temp_N[order(temp_N$gene_name, temp_N$gid, -temp_N$trx_len_N),]
#temp_H <- temp_H[order(temp_H$gene_name, temp_H$gid, -temp_H$trx_len_H),]
#temp_N_mono <- temp_N_mono[order(temp_N_mono$gene_name, temp_N_mono$gid, -temp_N_mono$trx_len_NHDF_DRACH_mono_dist),]
#temp_H_mono <- temp_H_mono[order(temp_H_mono$gene_name, temp_H_mono$gid, -temp_H_mono$trx_len_HD_DRACH_mono_dist),]

#temp_N <- temp_N[!duplicated(temp_N$gene_name),]
#temp_H <- temp_H[!duplicated(temp_H$gene_name),]
#temp_N_mono <- temp_N_mono[!duplicated(temp_N_mono$gene_name),]
#temp_H_mono <- temp_H_mono[!duplicated(temp_H_mono$gene_name),]


# limit m6a data to one transcript per gene (longest)
#NHDF_DRACH_dist <- NHDF_DRACH_dist[NHDF_DRACH_dist$refseqID %in% temp_N$gid,]
#HD_DRACH_dist <- HD_DRACH_dist[HD_DRACH_dist$refseqID %in% temp_H$gid,]
#NHDF_DRACH_mono_dist <- NHDF_DRACH_mono_dist[NHDF_DRACH_mono_dist$refseqID %in% temp_N_mono$gid,]
#HD_DRACH_mono_dist <- HD_DRACH_mono_dist[HD_DRACH_mono_dist$refseqID %in% temp_H_mono$gid,]



utr5.SF_N_D <- median(NHDF_DRACH_dist$utr5_size, na.rm = T)/median(NHDF_DRACH_dist$cds_size, na.rm = T)
utr3.SF_N_D <- median(NHDF_DRACH_dist$utr3_size, na.rm = T)/median(NHDF_DRACH_dist$cds_size, na.rm = T)

utr5.SF_H_D <- median(HD_DRACH_dist$utr5_size, na.rm = T)/median(HD_DRACH_dist$cds_size, na.rm = T)
utr3.SF_H_D <- median(HD_DRACH_dist$utr3_size, na.rm = T)/median(HD_DRACH_dist$cds_size, na.rm = T)

utr5.SF_N_D_mono <- median(NHDF_DRACH_mono_dist$utr5_size, na.rm = T)/median(NHDF_DRACH_mono_dist$cds_size, na.rm = T)
utr3.SF_N_D_mono <- median(NHDF_DRACH_mono_dist$utr3_size, na.rm = T)/median(NHDF_DRACH_mono_dist$cds_size, na.rm = T)

utr5.SF_H_D_mono <- median(HD_DRACH_mono_dist$utr5_size, na.rm = T)/median(HD_DRACH_mono_dist$cds_size, na.rm = T)
utr3.SF_H_D_mono <- median(HD_DRACH_mono_dist$utr3_size, na.rm = T)/median(HD_DRACH_mono_dist$cds_size, na.rm = T)

utr5.SF_N_a <- median(NHDF_a_dist$utr5_size, na.rm = T)/median(NHDF_a_dist$cds_size, na.rm = T)
utr3.SF_N_a <- median(NHDF_a_dist$utr3_size, na.rm = T)/median(NHDF_a_dist$cds_size, na.rm = T)

utr5.SF_H_a <- median(HD_a_dist$utr5_size, na.rm = T)/median(HD_a_dist$cds_size, na.rm = T)
utr3.SF_H_a <- median(HD_a_dist$utr3_size, na.rm = T)/median(HD_a_dist$cds_size, na.rm = T)

utr5.SF_N_a_mono <- median(NHDF_a_mono_dist$utr5_size, na.rm = T)/median(NHDF_a_mono_dist$cds_size, na.rm = T)
utr3.SF_N_a_mono <- median(NHDF_a_mono_dist$utr3_size, na.rm = T)/median(NHDF_a_mono_dist$cds_size, na.rm = T)

utr5.SF_H_a_mono <- median(HD_a_mono_dist$utr5_size, na.rm = T)/median(HD_a_mono_dist$cds_size, na.rm = T)
utr3.SF_H_a_mono <- median(HD_a_mono_dist$utr3_size, na.rm = T)/median(HD_a_mono_dist$cds_size, na.rm = T)


# assign the regions to new dataframes
utr5.NHDF_DRACH <- NHDF_DRACH_dist[NHDF_DRACH_dist$rel_location < 1, ]
cds.NHDF_DRACH <- NHDF_DRACH_dist [NHDF_DRACH_dist$rel_location < 2 & NHDF_DRACH_dist$rel_location >= 1, ]
utr3.NHDF_DRACH <- NHDF_DRACH_dist[NHDF_DRACH_dist$rel_location >= 2, ]

utr5.H_DRACH <- HD_DRACH_dist[HD_DRACH_dist$rel_location < 1, ]
cds.H_DRACH <- HD_DRACH_dist [HD_DRACH_dist$rel_location < 2 & HD_DRACH_dist$rel_location >= 1, ]
utr3.H_DRACH <- HD_DRACH_dist[HD_DRACH_dist$rel_location >= 2, ]

utr5.NHDF_DRACH_mono <- NHDF_DRACH_mono_dist[NHDF_DRACH_mono_dist$rel_location < 1, ]
cds.NHDF_DRACH_mono <- NHDF_DRACH_mono_dist [NHDF_DRACH_mono_dist$rel_location < 2 & NHDF_DRACH_mono_dist$rel_location >= 1, ]
utr3.NHDF_DRACH_mono <- NHDF_DRACH_mono_dist[NHDF_DRACH_mono_dist$rel_location >= 2, ]

utr5.HD_DRACH_mono <- HD_DRACH_mono_dist[HD_DRACH_mono_dist$rel_location < 1, ]
cds.HD_DRACH_mono <- HD_DRACH_mono_dist [HD_DRACH_mono_dist$rel_location < 2 & HD_DRACH_mono_dist$rel_location >= 1, ]
utr3.HD_DRACH_mono <- HD_DRACH_mono_dist[HD_DRACH_mono_dist$rel_location >= 2, ]

utr5.NHDF_a <- NHDF_a_dist[NHDF_a_dist$rel_location < 1, ]
cds.NHDF_a <- NHDF_a_dist [NHDF_a_dist$rel_location < 2 & NHDF_a_dist$rel_location >= 1, ]
utr3.NHDF_a <- NHDF_a_dist[NHDF_a_dist$rel_location >= 2, ]

utr5.H_a <- HD_a_dist[HD_a_dist$rel_location < 1, ]
cds.H_a <- HD_a_dist [HD_a_dist$rel_location < 2 & HD_a_dist$rel_location >= 1, ]
utr3.H_a <- HD_a_dist[HD_a_dist$rel_location >= 2, ]

utr5.NHDF_a_mono <- NHDF_a_mono_dist[NHDF_a_mono_dist$rel_location < 1, ]
cds.NHDF_a_mono <- NHDF_a_mono_dist [NHDF_a_mono_dist$rel_location < 2 & NHDF_a_mono_dist$rel_location >= 1, ]
utr3.NHDF_a_mono <- NHDF_a_mono_dist[NHDF_a_mono_dist$rel_location >= 2, ]

utr5.HD_a_mono <- HD_a_mono_dist[HD_a_mono_dist$rel_location < 1, ]
cds.HD_a_mono <- HD_a_mono_dist [HD_a_mono_dist$rel_location < 2 & HD_a_mono_dist$rel_location >= 1, ]
utr3.HD_a_mono <- HD_a_mono_dist[HD_a_mono_dist$rel_location >= 2, ]



#rescaling
utr5.NHDF_DRACH$rel_location <- rescale(utr5.NHDF_DRACH$rel_location, to = c(1-utr5.SF_N_D, 1), from = c(0,1))
utr3.NHDF_DRACH$rel_location <- rescale(utr3.NHDF_DRACH$rel_location, to = c(2, 2+utr3.SF_N_D), from = c(2,3))

utr5.H_DRACH$rel_location <- rescale(utr5.H_DRACH$rel_location, to = c(1-utr5.SF_H_D, 1), from = c(0,1))
utr3.H_DRACH$rel_location <- rescale(utr3.H_DRACH$rel_location, to = c(2, 2+utr3.SF_H_D), from = c(2,3))

utr5.NHDF_DRACH_mono$rel_location <- rescale(utr5.NHDF_DRACH_mono$rel_location, to = c(1-utr5.SF_N_D_mono, 1), from = c(0,1))
utr3.NHDF_DRACH_mono$rel_location <- rescale(utr3.NHDF_DRACH_mono$rel_location, to = c(2, 2+utr3.SF_N_D_mono), from = c(2,3))

utr5.HD_DRACH_mono$rel_location <- rescale(utr5.HD_DRACH_mono$rel_location, to = c(1-utr5.SF_H_D_mono, 1), from = c(0,1))
utr3.HD_DRACH_mono$rel_location <- rescale(utr3.HD_DRACH_mono$rel_location, to = c(2, 2+utr3.SF_H_D_mono), from = c(2,3))


utr5.NHDF_a$rel_location <- rescale(utr5.NHDF_a$rel_location, to = c(1-utr5.SF_N_a, 1), from = c(0,1))
utr3.NHDF_a$rel_location <- rescale(utr3.NHDF_a$rel_location, to = c(2, 2+utr3.SF_N_a), from = c(2,3))

utr5.H_a$rel_location <- rescale(utr5.H_a$rel_location, to = c(1-utr5.SF_H_a, 1), from = c(0,1))
utr3.H_a$rel_location <- rescale(utr3.H_a$rel_location, to = c(2, 2+utr3.SF_H_a), from = c(2,3))

utr5.NHDF_a_mono$rel_location <- rescale(utr5.NHDF_a_mono$rel_location, to = c(1-utr5.SF_N_a_mono, 1), from = c(0,1))
utr3.NHDF_a_mono$rel_location <- rescale(utr3.NHDF_a_mono$rel_location, to = c(2, 2+utr3.SF_N_a_mono), from = c(2,3))

utr5.HD_a_mono$rel_location <- rescale(utr5.HD_a_mono$rel_location, to = c(1-utr5.SF_H_a_mono, 1), from = c(0,1))
utr3.HD_a_mono$rel_location <- rescale(utr3.HD_a_mono$rel_location, to = c(2, 2+utr3.SF_H_a_mono), from = c(2,3))



N_D.metagene.coord <- c(utr5.NHDF_DRACH$rel_location, cds.NHDF_DRACH$rel_location, utr3.NHDF_DRACH$rel_location)
H_D.metagene.coord <- c(utr5.H_DRACH$rel_location, cds.H_DRACH$rel_location, utr3.H_DRACH$rel_location)
N_D_mono.metagene.coord <- c(utr5.NHDF_DRACH_mono$rel_location, cds.NHDF_DRACH_mono$rel_location, utr3.NHDF_DRACH_mono$rel_location)
H_D_mono.metagene.coord <- c(utr5.HD_DRACH_mono$rel_location, cds.HD_DRACH_mono$rel_location, utr3.HD_DRACH_mono$rel_location)

N_a.metagene.coord <- c(utr5.NHDF_a$rel_location, cds.NHDF_a$rel_location, utr3.NHDF_a$rel_location)
H_a.metagene.coord <- c(utr5.H_a$rel_location, cds.H_a$rel_location, utr3.H_a$rel_location)
N_a_mono.metagene.coord <- c(utr5.NHDF_a_mono$rel_location, cds.NHDF_a_mono$rel_location, utr3.NHDF_a_mono$rel_location)
H_a_mono.metagene.coord <- c(utr5.HD_a_mono$rel_location, cds.HD_a_mono$rel_location, utr3.HD_a_mono$rel_location)


#for plots
NHDF_all <- c(N_D.metagene.coord, N_a.metagene.coord)
Model <- c(rep("DRACH-context", length(N_D.metagene.coord)), 
           rep("All-context", length(N_a.metagene.coord))) 
df_NHDF <- data.frame(NHDF_all, Model)

NHDF_monoexonic <- c(N_D_mono.metagene.coord, N_a_mono.metagene.coord)
Model <- c(rep("DRACH-context", length(N_D_mono.metagene.coord)), 
           rep("All-context", length(N_a_mono.metagene.coord))) 
df_NHDF_mono <- data.frame(NHDF_monoexonic, Model)

HD10_6_all <- c(H_D.metagene.coord, H_a.metagene.coord)
Model <- c(rep("DRACH-context", length(H_D.metagene.coord)), 
           rep("All-context", length(H_a.metagene.coord))) 
df_HD <- data.frame(HD10_6_all, Model)

HD10_6_monoexonic <- c(H_D_mono.metagene.coord, H_a_mono.metagene.coord)
Model <- c(rep("DRACH-context", length(H_D_mono.metagene.coord)), 
           rep("All-context", length(H_a_mono.metagene.coord))) 
df_HD_mono <- data.frame(HD10_6_monoexonic, Model)


# Plots
NHDF_all_plot <- qplot(NHDF_all, data = df_NHDF, geom = "density", color = Model, size = I(1.5)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("All-context" = "darkred", "DRACH-context" = "red")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background shading
    axis.title = element_text(size = 21),  # Axis labels
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.text = element_text(size = 21),   # Axis tick labels
    axis.text.x = element_text(size = 25),   # X-axis tick labels
    axis.text.y = element_text(size = 25), 
    legend.text = element_text(size = 21), # Legend labels
    legend.title = element_text(size = 21) # Legend title
  )
NHDF_all_plot

NHDF_mono_plot <- qplot(NHDF_monoexonic, data = df_NHDF_mono, geom = "density", color = Model, size = I(1.5)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("All-context" = "darkred", "DRACH-context" = "red")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background shading
    axis.title = element_text(size = 21),  # Axis labels
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.text = element_text(size = 21),   # Axis tick labels
    axis.text.x = element_text(size = 25),   # X-axis tick labels
    axis.text.y = element_text(size = 25), 
    legend.text = element_text(size = 21), # Legend labels
    legend.title = element_text(size = 21) # Legend title
  )
NHDF_mono_plot


HD_all <- qplot(HD10_6_all, data = df_HD, geom = "density", color = Model, size = I(1.5)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("All-context" = "darkred", "DRACH-context" = "red")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background shading
    axis.title = element_text(size = 21),  # Axis labels
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.text = element_text(size = 21),   # Axis tick labels
    axis.text.x = element_text(size = 25),   # X-axis tick labels
    axis.text.y = element_text(size = 25), 
    legend.text = element_text(size = 21), # Legend labels
    legend.title = element_text(size = 21) # Legend title
  )
HD_all

HD_mono <- qplot(HD10_6_monoexonic, data = df_HD_mono, geom = "density", color = Model, size = I(1.5)) + 
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  geom_vline(xintercept = 1:2, col = "black") + 
  scale_color_manual(values = c("All-context" = "darkred", "DRACH-context" = "red")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove background shading
    axis.title = element_text(size = 21),  # Axis labels
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.text = element_text(size = 21),   # Axis tick labels
    axis.text.x = element_text(size = 25),   # X-axis tick labels
    axis.text.y = element_text(size = 25), 
    legend.text = element_text(size = 21), # Legend labels
    legend.title = element_text(size = 21) # Legend title
  )
HD_mono




combined_plot <- NHDF_all_plot + NHDF_mono_plot + HD_all + HD_mono + plot_layout(ncol = 2)
combined_plot




pdf(file = "NHDF_HD10.6_m6A_DRACH_and_allcontext_all_and_monoexonic_DMSO_v0.9.0_Filtered_stoichiometry_dist_hg38_10Stoich.pdf", width = 22, height = 15)
combined_plot
dev.off()



