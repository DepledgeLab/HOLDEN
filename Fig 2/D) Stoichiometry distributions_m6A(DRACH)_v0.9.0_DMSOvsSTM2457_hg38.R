library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

#setwd

#Filtering 20 reads
Unfiltered_reads_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif %>% filter(V11 > 0, V10 >=20)

Filtered_reads_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif.20reads_0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif %>% filter(V11 > 0, V10 >=20)


#filtering 10% stoichiometry
Unf_10_stoich_DMSO <- Unfiltered_reads_DMSO %>% filter(V11 >= 10)
Unf_10_stoich_STM <- Unfiltered_reads_STM %>% filter(V11 >= 10)
Filt_10_stoich_DMSO <- Filtered_reads_DMSO %>% filter(V11 >= 10)
Filt_10_stoich_STM <- Filtered_reads_STM %>% filter(V11 >= 10)



Unfiltered_10 <-ggplot() +
  geom_histogram(data = Unf_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unf_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Unfiltered, m6A >= 10%, DMSO(131,652) & STM2457(15,185)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, m6A >= 10%\nDMSO(131,652) & STM2457(15,185)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 13000)
Unfiltered_10
Filtered_10 <-ggplot() +
  geom_histogram(data = Filt_10_stoich_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filt_10_stoich_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Filtered, m6A >= 10%, DMSO(118,374) & STM2457(10,997)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, m6A >= 10%\nDMSO(118,374) & STM2457(10,997)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 13000)
Filtered_10


combined_plot <- Unfiltered_10 + Filtered_10  + plot_layout(ncol = 2)
combined_plot
pdf(file = "D) m6A_DRACH_v0.9.0_Unfiltered_Filtered_DMSOvsSTM2457_stoichiometry_10_dist_hg38.pdf", width = 15, height = 7.5)
combined_plot
dev.off()



