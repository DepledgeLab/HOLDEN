library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO <- pseU.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_STM <- pseU.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.noFilt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)

Filtered_reads_DMSO <- pseU.HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM <- pseU.HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.pseU.Filt.motif_20read_0stoich %>% filter(V11 > 0, V10 >=20)


#filtering 10% stoichiometry
Unf_10_stoich_DMSO <- Unfiltered_reads_DMSO %>% filter(V11 >= 10)
Unf_10_stoich_STM <- Unfiltered_reads_STM %>% filter(V11 >= 10)
Filt_10_stoich_DMSO <- Filtered_reads_DMSO %>% filter(V11 >= 10)
Filt_10_stoich_STM <- Filtered_reads_STM %>% filter(V11 >= 10)



#Plotting
Unfiltered <-ggplot() +
  geom_histogram(data = Unfiltered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Unfiltered,DMSO(3,435,705) & STM2457(5,491,596)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered\nDMSO(3,435,705) & STM2457(5,491,596)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3100000)
Unfiltered
Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Filtered, DMSO(564,062) & STM2457(1,004,940)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(564,062) & STM2457(1,004,940)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3100000)
Filtered
Zoomed_Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Filtered, DMSO(564,062) & STM2457(1,004,940)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(564,062) & STM2457(1,004,940)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 900000)
Zoomed_Filtered



Unfiltered_10 <-ggplot() +
  geom_histogram(data = Unf_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) + 
  geom_histogram(data = Unf_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Unfiltered, pseU >= 10%, DMSO(323,944) & STM2457(445,200)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, pseU >= 10%\nDMSO(323,944) & STM2457(445,200)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 120000)
Unfiltered_10
Filtered_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Filtered, pseU >= 10%, DMSO(7,035) & STM2457(10,397)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, pseU >= 10%\nDMSO(7,035) & STM2457(10,397)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 120000)
Filtered_10
Zoomed_Filtered_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#c7d6fc", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#0933a3", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "pseU stoichiometry", y = "Count", title = "Filtered, pseU >= 10%, DMSO(7,035) & STM2457(10,397)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, pseU >= 10%\nDMSO(7,035) & STM2457(10,397)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 3000)
Zoomed_Filtered_10

combined_plot <- Unfiltered + Filtered + Zoomed_Filtered + Unfiltered_10 + Filtered_10  + Zoomed_Filtered_10 + plot_layout(ncol = 3)
combined_plot
pdf(file = "A) HD_10_6_pseU_v0.9.0_Unfiltered_Filtered_DMSOvsSTM2457_stoichiometry_dist_alland10_hg38_supplementary.pdf", width = 26, height = 17)
combined_plot
dev.off()




