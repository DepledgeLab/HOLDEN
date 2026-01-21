library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif.stoich0  %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_STM <-  a.NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif.20reads.0stoich %>% filter(V11 > 0, V10 >=20)

Filtered_reads_DMSO <- a.NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20reads_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM <- a.NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20reads.stoich0 %>% filter(V11 > 0, V10 >=20)



#Plotting
Unfiltered <-ggplot() +
  geom_histogram(data = Unfiltered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all)  stoichiometry", y = "Count", title = "Unfiltered,DMSO(6,212,071) & STM2457(5,971,757)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered\nDMSO(6,212,071) & STM2457(5,971,757)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 4000000)
Unfiltered
Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, DMSO(1,875,411) & STM2457(1,529,028)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(1,875,411) & STM2457(1,529,028)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 4000000)
Filtered
Zoomed_Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A (all) stoichiometry", y = "Count", title = "Filtered, DMSO(1,875,411) & STM2457(1,529,028)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(1,875,411) & STM2457(1,529,028)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 1500000)
Zoomed_Filtered



combined_plot <- Unfiltered + Filtered + Zoomed_Filtered + plot_layout(ncol = 3)
combined_plot
pdf(file = "A) NHDF_m6A_all_context_v0.9.0_Unfiltered_Filtered_DMSOvsSTM2457_stoichiometry_dist_hg38_supplementary.pdf", width = 22, height = 7.5)
combined_plot
dev.off()


