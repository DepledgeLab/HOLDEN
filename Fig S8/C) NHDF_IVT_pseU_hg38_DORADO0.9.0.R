library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)


#setwd


#Filtering 20 reads all
Unfiltered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.noFilt.motif.20reads0stoich %>% filter (V4 %in% c("17802"))
Filtered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.Filt.motif.20reads0stoich %>% filter (V4 %in% c("17802"))

###filtering 10% stoichiometry-supplementary
Unfiltered_reads_IVT_all_10 <- Unfiltered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)
Filtered_reads_IVT_all_10 <- Filtered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)

Unfiltered_all <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, IVT(143,220 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, IVT(143,220 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 75000)
Unfiltered_all
Filtered_all <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU  stoichiometry", y = "Count", title = "Filtered, v0.9.0, IVT(24,573 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, IVT(24,573 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 75000)
Filtered_all


Unfiltered_10 <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, pseU >= 10%, IVT(14,315 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, pseU >= 10%\nIVT(14,315 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000)
Unfiltered_10
Filtered_10 <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU  stoichiometry", y = "Count", title = "Filtered, v0.9.0, pseU >= 10%, IVT(159 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, pseU >= 10%\nIVT(159 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000)
Filtered_10

combined_plot_pseu_10 <- Unfiltered_all + Filtered_all + Unfiltered_10 + Filtered_10 + plot_layout(ncol = 2)
combined_plot_pseu_10
pdf(file = "C) NHDF_IVT_pseU_v0.9.0_Unfiltered_Filtered_stoichiometry_dist_alland10_hg38.pdf", width = 17, height = 16)
combined_plot_pseu_10
dev.off()

