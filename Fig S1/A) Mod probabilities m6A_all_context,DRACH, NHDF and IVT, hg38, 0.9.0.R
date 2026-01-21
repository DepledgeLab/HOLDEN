#DMSO and STM - 0.9.0 -  all context
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-")& primary_base == "A")
A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-")& primary_base == "A")
A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-")& primary_base == "A")


A_merged_df_9 <- A_hg38_filtered_9_DMSO %>%
  inner_join(A_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_IVT_filtered_9, by = c("range_start", "range_end"))

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )

A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(A_plot_long_9$Sample)

custom_colors_A <- c("black", "black","grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_A) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "A(all) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 37, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
A_9


#DRACH

D_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-"))
D_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-"))
D_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-"))



D_merged_df_9 <- D_hg38_filtered_9_DMSO %>%
  inner_join(D_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(D_IVT_filtered_9, by = c("range_start", "range_end")) 

D_merged_df_9 <- D_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )


D_merged_df_9 <- D_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)


# Reshape the data to long format for ggplot2
D_plot_long_9 <- melt(D_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(D_plot_long_9$Sample)


custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
D_9 <- ggplot(D_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_A) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "A(DRACH) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 37, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
D_9

#combine the plots
combined_plot_9 <- (A_9 | D_9 | plot_spacer())
combined_plot_9

pdf(file = "A) A(all)_and_A(DRACH)in DMSO&STM in hg38&IVT,0.9.0,Profile Plot.pdf", width = 31, height = 9)
combined_plot_9
dev.off()

