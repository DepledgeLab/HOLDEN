library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO_7 <- NHDF_DMSO_48h_1.sup.m6A.trimAdapters.dorado.0.7.0.hg38.sorted.m6A.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_DMSO_8 <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.8.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_DMSO <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif.stoich0  %>% filter(V11 > 0, V10 >=20)


Filtered_reads_DMSO_7 <- NHDF_DMSO_48h_1.sup.m6A.trimAdapters.dorado.0.7.0.hg38.sorted.m6A.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_8 <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.8.0.hg38.sorted.m6A.inosine.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20reads_0stoich %>% filter(V11 > 0, V10 >=20)



#filtering 10% stoichiometry
Unf_10_stoich_DMSO_7 <- Unfiltered_reads_DMSO_7 %>% filter(V11 >= 10)
Unf_10_stoich_DMSO_8 <- Unfiltered_reads_DMSO_8 %>% filter(V11 >= 10)
Unf_10_stoich_DMSO_9 <- Unfiltered_reads_DMSO %>% filter(V11 >= 10)


Filt_10_stoich_DMSO_7 <- Filtered_reads_DMSO_7 %>% filter(V11 >= 10)
Filt_10_stoich_DMSO_8 <- Filtered_reads_DMSO_8 %>% filter(V11 >= 10)
Filt_10_stoich_DMSO_9 <- Filtered_reads_DMSO %>% filter(V11 >= 10)


#Plotting
Unfiltered <- ggplot() +
  geom_freqpoly(data = Unfiltered_reads_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unfiltered_reads_DMSO, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Unfiltered, v0.7.0, V0.8.0, v0.9.0") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 4000000)
Unfiltered

Filtered <- ggplot() +
  geom_freqpoly(data = Filtered_reads_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Filtered, v0.7.0, V0.8.0, v0.9.0") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 4000000)
Filtered


Zoomed_Filtered <- ggplot() +
  geom_freqpoly(data = Filtered_reads_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filtered_reads_DMSO, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Filtered, v0.7.0, V0.8.0, v0.9.0") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 1500000)
Zoomed_Filtered


Unfiltered_10 <- ggplot() +
  geom_freqpoly(data = Unf_10_stoich_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unf_10_stoich_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Unf_10_stoich_DMSO_9, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Unfiltered, v0.7.0, V0.8.0, v0.9.0\nm6A >= 10%") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 200000)
Unfiltered_10


Filtered_10 <- ggplot() +
  geom_freqpoly(data = Filt_10_stoich_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filt_10_stoich_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filt_10_stoich_DMSO_9, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Filtered, v0.7.0, V0.8.0, v0.9.0\nm6A >= 10%") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 200000)
Filtered_10


Zoomed_Filtered_10 <- ggplot() +
  geom_freqpoly(data = Filt_10_stoich_DMSO_7, aes(x = V11, color = "v0.7.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filt_10_stoich_DMSO_8, aes(x = V11, color = "v0.8.0"),
                binwidth = 2, size = 1.2) +
  geom_freqpoly(data = Filt_10_stoich_DMSO_9, aes(x = V11, color = "v0.9.0"),
                binwidth = 2, size = 1.2) +
  scale_color_manual(
    values = c(
      "v0.7.0" = "#FB8075",  # orange
      "v0.8.0" = "red",  # blue
      "v0.9.0" = "#B30000"   # green
    )
  ) +
  scale_x_continuous(limits = c(10, 100), breaks = seq(10, 100, 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count",
       title = "Filtered, v0.7.0, V0.8.0, v0.9.0\nm6A >= 10%") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),
    axis.text = element_text(size = 26, colour = "black"),
    axis.title = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  ) +
  ylim(0, 20000)
Zoomed_Filtered_10






combined_plot <- Unfiltered + Filtered + Zoomed_Filtered + Unfiltered_10 + Filtered_10  + Zoomed_Filtered_10 + plot_layout(ncol = 3)
combined_plot
pdf(file = "A) NHDF_m6A_all_context_v0.7.0_0.8.0_0.9.0_Unfiltered_Filtered_DMSO_stoichiometry_alland10_dist_hg38_supplementary.pdf", width = 27, height = 18)
combined_plot
dev.off()



