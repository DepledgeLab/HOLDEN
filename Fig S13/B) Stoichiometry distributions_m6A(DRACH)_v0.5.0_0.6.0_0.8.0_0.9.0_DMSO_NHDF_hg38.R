library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO_5 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.5.0.hg38.sorted.m6A.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_DMSO_6 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.6.0.hg38.sorted.m6A.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_DMSO_8 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.8.0.hg38.sorted.m6A.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_DMSO_9 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif %>% filter(V11 > 0, V10 >=20)

Filtered_reads_DMSO_5 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.5.0.hg38.sorted.m6A.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_6 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.6.0.hg38.sorted.m6A.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_8 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.8.0.hg38.sorted.m6A.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_9 <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif.20reads_0stoich %>% filter(V11 > 0, V10 >=20)



###filtering 10% stoichiometry-supplementary
Unfiltered_reads_DMSO_5_10 <- Unfiltered_reads_DMSO_5 %>% filter(V11 >= 10)
Unfiltered_reads_DMSO_6_10 <- Unfiltered_reads_DMSO_6 %>% filter(V11 >= 10)
Unfiltered_reads_DMSO_8_10 <- Unfiltered_reads_DMSO_8 %>% filter(V11 >= 10)
Unfiltered_reads_DMSO_9_10 <- Unfiltered_reads_DMSO_9 %>% filter(V11 >= 10)

Filtered_reads_DMSO_5_10 <- Filtered_reads_DMSO_5 %>% filter(V11 >= 10)
Filtered_reads_DMSO_6_10 <- Filtered_reads_DMSO_6 %>% filter(V11 >= 10)
Filtered_reads_DMSO_8_10 <- Filtered_reads_DMSO_8 %>% filter(V11 >= 10)
Filtered_reads_DMSO_9_10 <- Filtered_reads_DMSO_9 %>% filter(V11 >= 10)



#Plotting
Unfiltered <- ggplot() +
  geom_freqpoly(data = Unfiltered_reads_DMSO_5, aes(x = V11, color = "v0.5.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_6, aes(x = V11, color = "v0.6.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_9, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.5.0" = "#FDBCA0",
      "v0.6.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "m6A (DRACH) stoichiometry", y = "Count",
       title = "Unfiltered, v0.5.0, V0.6.0, v0.8.0, v0.9.0") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 250000)
Unfiltered

Filtered <- ggplot() +
  geom_freqpoly(data = Filtered_reads_DMSO_5, aes(x = V11, color = "v0.5.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_6, aes(x = V11, color = "v0.6.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_9, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.5.0" = "#FDBCA0",
      "v0.6.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "m6A (DRACH) stoichiometry", y = "Count",
       title = "Filtered, v0.5.0, V0.6.0, v0.8.0, v0.9.0") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 250000)
Filtered


Unfiltered_10 <- ggplot() +
  geom_freqpoly(data = Unfiltered_reads_DMSO_5_10, aes(x = V11, color = "v0.5.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_6_10, aes(x = V11, color = "v0.6.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_8_10, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_9_10, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.5.0" = "#FDBCA0",
      "v0.6.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  labs(x = "m6A (DRACH) stoichiometry", y = "Count",
       title = "Unfiltered, v0.5.0, V0.6.0, v0.8.0, v0.9.0\nm6A >= 10%") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 15000)
Unfiltered_10


Filtered_10 <- ggplot() +
  geom_freqpoly(data = Filtered_reads_DMSO_5_10, aes(x = V11, color = "v0.5.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_6_10, aes(x = V11, color = "v0.6.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_8_10, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_9_10, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.5.0" = "#FDBCA0",
      "v0.6.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  labs(x = "m6A (DRACH) stoichiometry", y = "Count",
       title = "Filtered, v0.5.0, V0.6.0, v0.8.0, v0.9.0\nm6A >= 10%") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 15000)
Filtered_10





combined_plot <- Unfiltered + Filtered + Unfiltered_10 + Filtered_10 + plot_layout(ncol = 2)
combined_plot
pdf(file = "B) NHDF_m6A_DRACH_v0.5.0_0.6.0_0.8.0_0.9.0_Unfiltered_Filtered_DMSO_stoichiometry_alland10_dist_hg38_supplementary.pdf", width = 17, height = 17)
combined_plot
dev.off()



