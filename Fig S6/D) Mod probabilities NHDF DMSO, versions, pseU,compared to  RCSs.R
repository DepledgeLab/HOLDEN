#DMSO  - pseU - 7,8,9
#Generating the merged dataframe with ranges and fractions
library(dplyr)
pseU_hg38_7 <- NHDF_DMSO_48h_1.sup.pseU.trimAdapters.dorado.0.7.0_hg38_probabilities %>% filter (code %in% c("17802"))
pseU_rcs_7 <- NHDF_DMSO_48h_1.sup.pseU.trimAdapters.dorado.0.7.0_rcs_probabilities %>% filter (code %in% c("17802"))
pseU_hg38_8 <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.8.0_hg38_probabilities %>% filter (code %in% c("17802")& primary_base == "T")
pseU_rcs_8 <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.8.0_rcs_probabilities %>% filter (code %in% c("17802")& primary_base == "T")
pseU_hg38_9 <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("17802")& primary_base == "T")
pseU_rcs_9 <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("17802")& primary_base == "T")


pseU_merged_df_all <- pseU_hg38_7 %>%
  inner_join(pseU_hg38_8, by = c("range_start", "range_end")) %>%
  inner_join(pseU_hg38_9, by = c("range_start", "range_end")) %>%
  inner_join(pseU_rcs_7, by = c("range_start", "range_end")) %>%
  inner_join(pseU_rcs_8, by = c("range_start", "range_end")) %>%
  inner_join(pseU_rcs_9, by = c("range_start", "range_end"))


pseU_merged_df_all <- pseU_merged_df_all %>%
  rename(
    NHDF_0.7.0 = frac.x,
    NHDF_0.8.0 = frac.y,
    NHDF_0.9.0 = frac.x.x,
    RCS_0.7.0 = frac.y.y,
    RCS_0.8.0 = frac.x.x.x,
    RCS_0.9.0 = frac.y.y.y
  )


pseU_merged_df_all <- pseU_merged_df_all %>% select(range_start, NHDF_0.7.0, NHDF_0.8.0, NHDF_0.9.0, RCS_0.7.0, RCS_0.8.0, RCS_0.9.0)

write.table(pseU_merged_df_all, "NHDF_DMSO_RCS_pseU_hg38&rcs_0.7_0.8_0.9_merged_df_alternative.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
pseU_merged_long_all <- melt(pseU_merged_df_all, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(pseU_merged_long_all$Sample)

custom_colors_pseU <- c('#a0c3ff','#669fff', "#346ccc", 'grey','grey52', "grey20")

custom_linetypes <- c("solid", "solid", "solid", "solid", "solid", "solid")

# Create the plots
pseU_all <-ggplot(pseU_merged_long_all, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_pseU) +  # Apply custom colors
  labs(x = "modification probability", y = "pseU fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.015)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
pseU_all
#pseU from 0.9 to 1.0 probability
pseU_merged_long_all_filt <- pseU_merged_long_all %>% filter(range_start >= 0.95)

pseU_all_filtered <-ggplot(pseU_merged_long_all_filt, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_pseU) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "modification probability", y = "pseU fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.01)) +  # Set limits and format y-axis
  scale_x_continuous(breaks = seq(0.95, 1, by = 0.01), expand = c(0.003, 0.004)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
pseU_all_filtered
#combine the plots
combined_plot_u <- ( pseU_all | pseU_all_filtered)
combined_plot_u

pdf(file = "C) pseU in DMSO in hg38&rcs,7,8,9 versions,Profile Plot.pdf", width = 14, height = 5)
combined_plot_u
dev.off()
