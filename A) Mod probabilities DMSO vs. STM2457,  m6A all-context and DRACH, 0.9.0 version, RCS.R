#DMSO and STM - 0.9.0 - m6A all-context
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a")& primary_base == "A")
A_rcs_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_rcs_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("a")& primary_base == "A")
A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a")& primary_base == "A")
A_rcs_filtered_9_STM <- NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_rcs_filtered_9_STM <- NHDF_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("a")& primary_base == "A")


A_merged_df_9 <- A_hg38_filtered_9_DMSO %>%
  inner_join(A_rcs_filtered_9_DMSO, by = c("range_start", "range_end")) %>%
  inner_join(A_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_rcs_filtered_9_STM, by = c("range_start", "range_end"))

m6A_merged_df_9 <- m6A_hg38_filtered_9_DMSO %>%
  inner_join(m6A_rcs_filtered_9_DMSO, by = c("range_start", "range_end")) %>%
  inner_join(m6A_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(m6A_rcs_filtered_9_STM, by = c("range_start", "range_end"))

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    DMSO_RCS = frac.y,
    STM2457 = frac.x.x,
    STM2457_RCS = frac.y.y,
  )
m6A_merged_df_9 <- m6A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    DMSO_RCS = frac.y,
    STM2457 = frac.x.x,
    STM2457_RCS = frac.y.y,
  )
A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, DMSO_RCS, STM2457, STM2457_RCS)
m6A_merged_df_9 <- m6A_merged_df_9 %>% select(range_start, DMSO, DMSO_RCS, STM2457, STM2457_RCS)

write.table(A_merged_df_9, "NHDF_DMSO_STM_A_hg38&rcs_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(m6A_merged_df_9, "NHDF_DMSO_STM_m6A_hg38&rcs_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Series", value.name = "Value")
m6A_plot_long_9 <- melt(m6A_merged_df_9, id.vars = "range_start", variable.name = "Series", value.name = "Value")

unique(A_plot_long_9$Series)
unique(m6A_plot_long_9$Series)

custom_colors_A <- c("black", "grey52","black", "grey52")
custom_colors_m6A <- c("red", "grey52", "red", "grey52")

custom_linetypes <- c("solid", "solid", "dashed", "dashed")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_A) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
A_9
m6A_9 <-ggplot(m6A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_m6A) +  # Apply custom colors
  labs(x = "modification probability", y = "m6A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.62)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m6A_9
#m6A from 0.9 to 1.0 probability
m6A_plot_long_9_filtered <- m6A_plot_long_9 %>% filter(range_start >= 0.95)

m6A_9_filtered <-ggplot(m6A_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_m6A) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "modification probability", y = "m6A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  # Set limits and format y-axis
  scale_x_continuous(breaks = seq(0.95, 1, by = 0.01), expand = c(0.003, 0.004)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m6A_9_filtered
#combine the plots
combined_plot_9 <- (A_9 | m6A_9 | m6A_9_filtered)
combined_plot_9

pdf(file = "D) NHDF,A&m6A all-context in DMSO&STM in hg38&rcs,0.9.0,Profile Plot.pdf", width = 22, height = 5)
combined_plot_9
dev.off()



#DMSO and STM - 0.9.0 - m6A Drach
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-"))
m6A_hg38_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a"))
A_rcs_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("-"))
m6A_rcs_filtered_9_DMSO <- NHDF_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("a"))
A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("-"))
m6A_hg38_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_hg38_probabilities %>% filter (code %in% c("a"))
A_rcs_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("-"))
m6A_rcs_filtered_9_STM <- NHDF_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_rcs_probabilities %>% filter (code %in% c("a"))



A_merged_df_9 <- A_hg38_filtered_9_DMSO %>%
  inner_join(A_rcs_filtered_9_DMSO, by = c("range_start", "range_end")) %>%
  inner_join(A_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_rcs_filtered_9_STM, by = c("range_start", "range_end"))

m6A_merged_df_9 <- m6A_hg38_filtered_9_DMSO %>%
  inner_join(m6A_rcs_filtered_9_DMSO, by = c("range_start", "range_end")) %>%
  inner_join(m6A_hg38_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(m6A_rcs_filtered_9_STM, by = c("range_start", "range_end"))

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    DMSO_RCS = frac.y,
    STM2457 = frac.x.x,
    STM2457_RCS = frac.y.y,
  )
m6A_merged_df_9 <- m6A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    DMSO_RCS = frac.y,
    STM2457 = frac.x.x,
    STM2457_RCS = frac.y.y,
  )
A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, DMSO_RCS, STM2457, STM2457_RCS)
m6A_merged_df_9 <- m6A_merged_df_9 %>% select(range_start, DMSO, DMSO_RCS, STM2457, STM2457_RCS)

write.table(A_merged_df_9, "NHDF_DMSO_STM_A_DRACH_hg38&rcs_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(m6A_merged_df_9, "NHDF_DMSO_STM_m6A_DRACH_hg38&rcs_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Series", value.name = "Value")
m6A_plot_long_9 <- melt(m6A_merged_df_9, id.vars = "range_start", variable.name = "Series", value.name = "Value")

unique(A_plot_long_9$Series)
unique(m6A_plot_long_9$Series)

custom_colors_A <- c("black", "grey52","black", "grey52")
custom_colors_m6A_DRACH <- c("red", "grey52", "red", "grey52")

custom_linetypes <- c("solid", "solid", "dashed", "dashed")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_A) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
A_9
m6A_9 <-ggplot(m6A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_m6A_DRACH) +  # Apply custom colors
  labs(x = "modification probability", y = "m6A(DRACH) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.62)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m6A_9
#m6A from 0.9 to 1.0 probability
m6A_plot_long_9_filtered <- m6A_plot_long_9 %>% filter(range_start >= 0.95)

m6A_9_filtered <-ggplot(m6A_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Series", group = "Series", linetype = "Series")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_m6A_DRACH) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "modification probability", y = "m6A(DRACH) fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.62)) +  # Set limits and format y-axis
  scale_x_continuous(breaks = seq(0.95, 1, by = 0.01), expand = c(0.003, 0.004)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m6A_9_filtered
#combine the plots
combined_plot_9 <- (A_9 | m6A_9 | m6A_9_filtered)
combined_plot_9

pdf(file = "D) NHDF, A&m6A(DRACH) in DMSO&STM in hg38&rcs,0.9.0,Profile Plot.pdf", width = 22, height = 5)
combined_plot_9
dev.off()



