library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)


#setwd

#All context m6A
#Filtering 20 reads all
Unfiltered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A_inosine.noFilt.motif.20reads0stoich %>% filter(V11 > 0, V10 >=20)
Filtered_reads_IVT_all <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A_inosine.Filt.motif.20reads0stoich %>% filter(V11 > 0, V10 >=20)
Unfiltered_reads_IVT_all <- Unfiltered_reads_IVT_all %>% filter (V4 %in% c("a"))
Filtered_reads_IVT_all <- Filtered_reads_IVT_all %>% filter (V4 %in% c("a"))

###filtering 10% stoichiometry-supplementary
Unfiltered_reads_IVT_all_10 <- Unfiltered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)
Filtered_reads_IVT_all_10 <- Filtered_reads_IVT_all %>% filter(V11 > 10, V10 >=20)



Unfiltered_all <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, IVT(163,150 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, IVT(163,150 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 75000)
Unfiltered_all
Filtered_all <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A  stoichiometry", y = "Count", title = "Filtered, v0.9.0, IVT(32,876 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, IVT(32,876 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 75000)
Filtered_all


Unfiltered_all_10 <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, m6A >= 10%, IVT(14,145 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, m6A >= 10%\nIVT(14,145 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500)
Unfiltered_all_10
Filtered_all_10 <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_all_10, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A  stoichiometry", y = "Count", title = "Filtered, v0.9.0, m6A >= 10%, IVT(264 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, m6A >= 10%\nIVT(264 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 2500)
Filtered_all_10


combined_plot_all <- Unfiltered_all + Filtered_all +  Unfiltered_all_10 + Filtered_all_10 + plot_layout(ncol = 2)
combined_plot_all 
pdf(file = "A) m6A_allcontext_0.9.0_Unfiltered_Filtered_NHDF_IVT_stoichiometry _alland10_dist_hg38_supplementary.pdf", width = 17, height = 16)
combined_plot_all
dev.off()



#DRACH context m6A
#Filtering 20 reads DRACH
Unfiltered_reads_IVT_DRACH <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif %>% filter(V11 > 0, V10 >=20)
Filtered_reads_IVT_DRACH <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif %>% filter(V11 > 0, V10 >=20)

###filtering 10% stoichiometry-supplementary
Unfiltered_reads_IVT_DRACH_10 <- Unfiltered_reads_IVT_DRACH %>% filter(V11 > 10, V10 >=20)
Filtered_reads_IVT_DRACH_10 <- Filtered_reads_IVT_DRACH %>% filter(V11 > 10, V10 >=20)



#Plotting
Unfiltered_DRACH <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_DRACH, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, IVT(9,335 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, IVT(9,335 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 6000)
Unfiltered_DRACH
Filtered_DRACH <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_DRACH, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Filtered, v0.9.0, IVT(2,195 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, IVT(2,195 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 6000)
Filtered_DRACH

###Plotting for supplementary figures 
Unfiltered_DRACH_10 <-ggplot() +
  geom_histogram(data = Unfiltered_reads_IVT_DRACH_10, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Unfiltered, v0.9.0, m6A >= 10%, IVT(93 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Unfiltered, v0.9.0, m6A >= 10%\nIVT(93 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30)
Unfiltered_DRACH_10
Filtered_DRACH_10 <-ggplot() +
  geom_histogram(data = Filtered_reads_IVT_DRACH_10, aes(x = V11), 
                 binwidth = 2, fill = "#f3511f", alpha = 0.9, color = "black", boundary = 0) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(x = "m6A(DRACH)  stoichiometry", y = "Count", title = "Filtered, v0.9.0, m6A >= 10%, IVT(3 sites)") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  ggtitle("Filtered, v0.9.0, m6A >= 10%\nIVT(3 sites)") +
  theme(
    plot.title = element_text(size = 27, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 29),  # Axis text size
    axis.title = element_text(size = 29)  # Axis title size
  ) +
  ylim(0, 30)
Filtered_DRACH_10


combined_plot_DRACH <- Unfiltered_DRACH + Filtered_DRACH + Unfiltered_DRACH_10 + Filtered_DRACH_10 + plot_layout(ncol = 2)
combined_plot_DRACH 
pdf(file = "B) m6A_DRACH_0.9.0_Unfiltered_Filtered_NHDF_IVT_stoichiometry_dist_alland10_hg38_supplementary.pdf", width = 17, height = 16)
combined_plot_DRACH
dev.off()



