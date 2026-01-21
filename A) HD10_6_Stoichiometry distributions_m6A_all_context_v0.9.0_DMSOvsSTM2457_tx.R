library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO <- a.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.inosine.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_STM <-  a.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.inosine.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

Filtered_reads_DMSO <- a.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.inosine.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM <- a.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.gencode_v47.sorted.m6A.inosine.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)


#filtering 10% stoichiometry
Unf_10_stoich_DMSO <- Unfiltered_reads_DMSO %>% filter(V11 >= 10)
Unf_10_stoich_STM <- Unfiltered_reads_STM %>% filter(V11 >= 10)
Filt_10_stoich_DMSO <- Filtered_reads_DMSO %>% filter(V11 >= 10)
Filt_10_stoich_STM <- Filtered_reads_STM %>% filter(V11 >= 10)



#Plotting
Unfiltered <-ggplot() +
  geom_histogram(data = Unfiltered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Unfiltered,DMSO(3,918,398) & STM2457(6,029,384)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered\nDMSO(3,918,398) & STM2457(6,029,384)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000000)
Unfiltered
Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, DMSO(980,651) & STM2457(1,244,304)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(980,651) & STM2457(1,244,304)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000000)
Filtered
Zoomed_Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, DMSO(980,651) & STM2457(1,244,304)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(980,651) & STM2457(1,244,304)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 1000000)
Zoomed_Filtered



Unfiltered_10 <-ggplot() +
  geom_histogram(data = Unf_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unf_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Unfiltered, m6A >= 10%, DMSO(385,817) & STM2457(332,989)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, m6A >= 10%\nDMSO(385,817) & STM2457(332,989)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 120000)
Unfiltered_10
Filtered_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, m6A >= 10%, DMSO(128,959) & STM2457(19,119)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, m6A >= 10%\nDMSO(128,959) & STM2457(19,119)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 120000)
Filtered_10
Zoomed_Filtered_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, m6A >= 10%, DMSO(128,959) & STM2457(19,119)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, m6A >= 10%\nDMSO(128,959) & STM2457(19,119)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 15000)
Zoomed_Filtered_10

combined_plot <- Unfiltered + Filtered + Zoomed_Filtered + Unfiltered_10 + Filtered_10  + Zoomed_Filtered_10 + plot_layout(ncol = 3)
combined_plot
pdf(file = "A) HD_10_6_m6A_all_context_v0.9.0_Unfiltered_Filtered_DMSOvsSTM2457_stoichiometry_dist_alland10_tx_supplementary.pdf", width = 26, height = 17)
combined_plot
dev.off()


