library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd


#Inosine, 20 reads
Unfiltered_reads_DMSO_I <- Inosine.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Filtered_reads_DMSO_I <- Inosine.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine99.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

Unfiltered_reads_STM_I <- Inosine.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif_20read_0stoich  %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM_I <- Inosine.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine99.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

#filtering 10% stoichiometry
Unf_10_stoich_DMSO_I <- Unfiltered_reads_DMSO_I %>% filter(V11 >= 10)
Filt_10_stoich_DMSO_I <- Filtered_reads_DMSO_I %>% filter(V11 >= 10)
Unf_10_stoich_STM_I <- Unfiltered_reads_STM_I %>% filter(V11 >= 10)
Filt_10_stoich_STM_I <- Filtered_reads_STM_I %>% filter(V11 >= 10)

#Plotting
Unfiltered_reads_I <-ggplot() +
  geom_histogram(data = Unfiltered_reads_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Unfiltered, DMSO(2,159,279) vs. STM2457(3,474,528)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered\nDMSO(2,159,279) vs. STM2457(3,474,528)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500000)
Unfiltered_reads_I
Filtered_reads_I <-ggplot() +
  geom_histogram(data = Filtered_reads_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, DMSO(347,179) vs. STM2457(536,809)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(347,179) vs. STM2457(536,809)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500000)
Filtered_reads_I
Zoomed_Filtered_I <-ggplot() +
  geom_histogram(data = Filtered_reads_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, DMSO(347,179) vs. STM2457(536,809)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(347,179) vs. STM2457(536,809)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 500000)
Zoomed_Filtered_I


Unfiltered_reads_I_10 <-ggplot() +
  geom_histogram(data = Unf_10_stoich_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unf_10_stoich_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill = "purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Unfiltered, I >= 10%, DMSO(88,633) vs. STM2457(109,925)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, I >= 10%\nDMSO(88,633) vs. STM2457(109,925)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30000)
Unfiltered_reads_I_10
Filtered_reads_I_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, I >= 10%, DMSO(5,614) vs. STM2457(14,288)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, I >= 10%\nDMSO(5,614) vs. STM2457(14,288)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30000)
Filtered_reads_I_10
Zoomed_Filtered_I_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM_I, aes(x = V11), 
                 binwidth = 2, fill = "#d3c7fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO_I, aes(x = V11), 
                 binwidth = 2, fill ="purple", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "I stoichiometry", y = "Count", title = "Filtered, I >= 10%, DMSO(5,614) vs. STM2457(14,288)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, I >= 10%\nDMSO(5,614) vs. STM2457(14,288)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500)
Zoomed_Filtered_I_10


combined_plot <- Unfiltered_reads_I + Filtered_reads_I + Zoomed_Filtered_I + Unfiltered_reads_I_10 + Filtered_reads_I_10  + Zoomed_Filtered_I_10 + plot_layout(ncol = 3)
combined_plot
pdf(file = "B) HD10_6_I_v0.9.0_Unfiltered_Filtered_DMSOvsSTM_stoichiometry_dist_hg38_inc.10stoich_supplementary.pdf", width = 26, height = 17)
combined_plot
dev.off()



