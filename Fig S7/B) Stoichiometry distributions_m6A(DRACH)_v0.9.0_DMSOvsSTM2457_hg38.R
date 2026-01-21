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


#Plotting
Unfiltered <-ggplot() +
  geom_histogram(data = Unfiltered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Unfiltered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Unfiltered,DMSO(539,170) & STM2457(397,655)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered\nDMSO(539,170) & STM2457(397,655)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 300000)
Unfiltered
Filtered <-ggplot() +
  geom_histogram(data = Filtered_reads_DMSO, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  geom_histogram(data = Filtered_reads_STM, aes(x = V11), 
                 binwidth = 2, fill = "#fbc2b1", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Filtered, DMSO(390,729) & STM2457(178,666)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered\nDMSO(390,729) & STM2457(178,666)") +
  theme(
    plot.title = element_text(size = 22, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 26),  # Axis text size
    axis.title = element_text(size = 26)  # Axis title size
  ) +
  ylim(0, 300000)
Filtered



combined_plot <- Unfiltered + Filtered  + plot_layout(ncol = 2)
combined_plot
pdf(file = "B) m6A_DRACH_v0.9.0_Unfiltered_Filtered_DMSOvsSTM2457_stoichiometry_dist_hg38.pdf", width = 15, height = 7.5)
combined_plot
dev.off()




