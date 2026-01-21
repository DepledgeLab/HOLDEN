#DMSO and STM & IVT - 0.9.0
#Generating the merged dataframe with ranges and fractions
library(dplyr)
m6A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a"))
m6A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a"))
m6A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a"))


#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)


# Create the plots
m6A_IVT <- ggplot(m6A_IVT_filtered_9, aes(x = range_start, y = frac)) +
  geom_col(fill = "grey52", color = "grey52") +
  labs(x = "modification probability (IVT)", y = "m6A(all) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 0.62)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
m6A_IVT

m6A_DMSO <- ggplot(m6A_hg38_filtered_9_DMSO, aes(x = range_start, y = frac)) +
  geom_col(fill = "red", color = "red") +
  labs(x = "modification probability (DMSO)", y = "m6A(all) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 0.62)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
m6A_DMSO

m6A_STM <- ggplot(m6A_hg38_filtered_9_STM, aes(x = range_start, y = frac)) +
  geom_col(fill = "red4", color = "red4") +
  labs(x = "modification probability (STM2457)", y = "m6A(all) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 0.62)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
m6A_STM


#combine the plots
combined_plot_9 <- (m6A_IVT | m6A_DMSO | m6A_STM)
combined_plot_9

pdf(file = "B) m6A(DRACH) in DMSO&STM and IVT,0.9.0,Profile Plot.pdf", width = 20, height = 6)
combined_plot_9
dev.off()
