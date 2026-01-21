library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#9, 20 reads
Unfiltered_reads_DMSO_9 <- Inosine.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_9 <- Inosine.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine99.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

Unfiltered_reads_STM_9 <- Inosine.NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM_9 <- Inosine.NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine99.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

#filtering 10% stoichiometry
Unf_10_stoich_DMSO_9 <- Unfiltered_reads_DMSO_9 %>% filter(V11 >= 10)
Filt_10_stoich_DMSO_9 <- Filtered_reads_DMSO_9 %>% filter(V11 >= 10)
Unf_10_stoich_STM_9 <- Unfiltered_reads_STM_9 %>% filter(V11 >= 10)
Filt_10_stoich_STM_9 <- Filtered_reads_STM_9 %>% filter(V11 >= 10)

#Plotting
Unfiltered_reads_9 <-ggplot() +
  geom_histogram(data = Unfiltered_reads_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, DMSO(3,667,939) vs. STM2457(3,686,141)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0\nDMSO(3,667,939) vs. STM2457(3,686,141)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000000)
Unfiltered_reads_9
Filtered_reads_9 <-ggplot() +
  geom_histogram(data = Filtered_reads_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, v0.9.0, DMSO(694,467) vs. STM2457(730,082)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0\nDMSO(694,467) vs. STM2457(730,082)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000000)
Filtered_reads_9
Zoomed_Filtered_9 <-ggplot() +
  geom_histogram(data = Filtered_reads_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, v0.9.0, DMSO(694,467) vs. STM2457(730,082)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0\nDMSO(694,467) vs. STM2457(730,082)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 700000)
Zoomed_Filtered_9


Unfiltered_reads_9_10 <-ggplot() +
  geom_histogram(data = Unf_10_stoich_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unf_10_stoich_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, I >= 10%, DMSO(116,548) vs. STM2457(119,832)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, I >= 10%\nDMSO(116,548) vs. STM2457(119,832)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30000)
Unfiltered_reads_9_10
Filtered_reads_9_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, v0.9.0, I >= 10%, DMSO(9,722) vs. STM2457(14,263)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, I >= 10%\nDMSO(9,722) vs. STM2457(14,263)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30000)
Filtered_reads_9_10
Zoomed_Filtered_9_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM_9, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO_9, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, v0.9.0, I >= 10%, DMSO(9,722) vs. STM2457(14,263)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, I >= 10%\nDMSO(9,722) vs. STM2457(14,263)") +
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500)
Zoomed_Filtered_9_10



combined_plot <- Unfiltered_reads_9 + Filtered_reads_9 + Zoomed_Filtered_9 + Unfiltered_reads_9_10 + Filtered_reads_9_10 + Zoomed_Filtered_9_10 + plot_layout(ncol = 3)
combined_plot
pdf(file = "B) NHDF_Inosine_v0.9.0_Unfiltered_Filtered_DMSOvsSTM_stoichiometry_dist_hg38_inc.10stoich_supplementary.pdf", width = 26, height = 17)
combined_plot
dev.off()

