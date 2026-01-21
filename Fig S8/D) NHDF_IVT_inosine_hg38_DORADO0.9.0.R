library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)


#setwd


#Filtering 20 reads all
Unfiltered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A_inosine.noFilt.motif.20reads0stoich %>% filter (V4 %in% c("17596"))
Filtered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A_inosine.Filt.motif.20reads0stoich %>% filter (V4 %in% c("17596"))

###filtering 10% stoichiometry-supplementary
Unfiltered_reads_IVT_all_10 <- Unfiltered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)
Filtered_reads_IVT_all_10 <- Filtered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)

Unfiltered_all <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Inosine  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, IVT(84,547 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, IVT(84,547 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 50000)
Unfiltered_all
Filtered_all <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Inosine  stoichiometry", y = "Count", title = "Filtered, v0.9.0, IVT(18,129 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, IVT(18,129 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 50000)
Filtered_all


Unfiltered_10 <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Inosine  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, Inosine >= 10%, IVT(9,244 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, Inosine >= 10%\nIVT(9,244 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2000)
Unfiltered_10
Filtered_10 <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "Inosine  stoichiometry", y = "Count", title = "Filtered, v0.9.0, Inosine >= 10%, IVT(150 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, Inosine >= 10%\nIVT(150 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2000)
Filtered_10

combined_plot_Inosine_10 <- Unfiltered_all + Filtered_all + Unfiltered_10  + Filtered_10 +  plot_layout(ncol = 2)
combined_plot_Inosine_10
pdf(file = "D) NHDF_IVT_Inosine_v0.9.0_Unfiltered_Filtered_stoichiometry_dist_alland10_hg38.pdf", width = 17, height = 16)
combined_plot_Inosine_10
dev.off()
