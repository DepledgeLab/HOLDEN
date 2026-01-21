#DMSO and STM - 0.9.0 - m6A all context
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a")& primary_base == "A")
A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a")& primary_base == "A")
A_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
m6A_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a")& primary_base == "A")


A_merged_df_9 <- A_tx_filtered_9_DMSO %>%
  inner_join(A_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_IVT_filtered_9, by = c("range_start", "range_end"))

m6A_merged_df_9 <- m6A_tx_filtered_9_DMSO %>%
  inner_join(m6A_IVT_filtered_9, by = c("range_start", "range_end")) %>%
  inner_join(m6A_tx_filtered_9_STM, by = c("range_start", "range_end"))

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
m6A_merged_df_9 <- m6A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )

A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)
m6A_merged_df_9 <- m6A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

write.table(A_merged_df_9, "HD10.6_DMSO_STM_A_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(m6A_merged_df_9, "HD10.6_DMSO_STM_m6A_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")
m6A_plot_long_9 <- melt(m6A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(A_plot_long_9$Sample)
unique(m6A_plot_long_9$Sample)

custom_colors_A <- c("black", "black","grey52")
custom_colors_m6A <- c("red", "red", "grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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
m6A_9 <-ggplot(m6A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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

m6A_9_filtered <-ggplot(m6A_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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

pdf(file = "A) HD10.6,A&m6A all-context in DMSO&STM in tx&IVT,0.9.0,Profile Plot.pdf", width = 20, height = 5)
combined_plot_9
dev.off()


#DMSO and STM & IVT - 0.9.0 - m6A DRACH
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-"))
m6A_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a"))
A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-"))
m6A_IVT_filtered_9 <- NHDF_polyA_TSO_IVT.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a"))
A_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-"))
m6A_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.m6A_DRACH.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("a"))


A_merged_df_9 <- A_tx_filtered_9_DMSO %>%
  inner_join(A_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_IVT_filtered_9, by = c("range_start", "range_end")) 

m6A_merged_df_9 <- m6A_tx_filtered_9_DMSO %>%
  inner_join(m6A_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(m6A_IVT_filtered_9, by = c("range_start", "range_end")) 

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
m6A_merged_df_9 <- m6A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )

A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)
m6A_merged_df_9 <- m6A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

write.table(A_merged_df_9, "HD10.6_DMSO_STM_A_DRACH_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(m6A_merged_df_9, "HD10.6_DMSO_STM_m6A_DRACH_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")
m6A_plot_long_9 <- melt(m6A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(A_plot_long_9$Sample)
unique(m6A_plot_long_9$Sample)

custom_colors_A <- c("black", "black","grey52")
custom_colors_m6A_DRACH <- c("red", "red", "grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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
m6A_9 <-ggplot(m6A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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

m6A_9_filtered <-ggplot(m6A_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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

pdf(file = "B) HD10.6,A&m6A(DRACH)in DMSO&STM in tx&IVT,0.9.0,Profile Plot.pdf", width = 20, height = 5)
combined_plot_9
dev.off()




#DMSO and STM - 0.9.0 - pseU
#Generating the merged dataframe with ranges and fractions
library(dplyr)
U_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "T")
pseU_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17802")& primary_base == "T")
U_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "T")
pseU_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17802")& primary_base == "T")
U_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "T")
pseU_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17802")& primary_base == "T")


U_merged_df_9 <- U_tx_filtered_9_DMSO %>%
  inner_join(U_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(U_IVT_filtered_9_DMSO, by = c("range_start", "range_end")) 

pseU_merged_df_9 <- pseU_tx_filtered_9_DMSO %>%
  inner_join(pseU_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(pseU_IVT_filtered_9_DMSO, by = c("range_start", "range_end"))

U_merged_df_9 <- U_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
pseU_merged_df_9 <- pseU_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
U_merged_df_9 <- U_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)
pseU_merged_df_9 <- pseU_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

write.table(U_merged_df_9, "HD10.6_DMSO_STM_U_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(pseU_merged_df_9, "HD10.6_DMSO_STM_pseU_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
U_plot_long_9 <- melt(U_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")
pseU_plot_long_9 <- melt(pseU_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")


custom_colors_U <- c("black", "black","grey52")
custom_colors_pseU <- c("blue", "blue", "grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
U_9 <- ggplot(U_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_U) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "U fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
U_9
pseU_9 <-ggplot(pseU_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_pseU) +  # Apply custom colors
  labs(x = "modification probability", y = "pseU fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.02)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
pseU_9
#pseU from 0.9 to 1.0 probability
pseU_plot_long_9_filtered <- pseU_plot_long_9 %>% filter(range_start >= 0.95)

pseU_9_filtered <-ggplot(pseU_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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
pseU_9_filtered
#combine the plots
combined_plot_9 <- (U_9 | pseU_9 | pseU_9_filtered)
combined_plot_9

pdf(file = "C) HD10.6,U&pseU in DMSO&STM in tx&IVT,0.9.0,Profile Plot.pdf", width = 20, height = 5)
combined_plot_9
dev.off()


#DMSO and STM - 0.9.0 - inosine
#Generating the merged dataframe with ranges and fractions
library(dplyr)
A_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
i_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17596")& primary_base == "A")
A_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
i_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17596")& primary_base == "A")
A_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "A")
i_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("17596")& primary_base == "A")

A_merged_df_9 <- A_tx_filtered_9_DMSO %>%
  inner_join(A_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(A_IVT_filtered_9_DMSO, by = c("range_start", "range_end"))

i_merged_df_9 <- i_tx_filtered_9_DMSO %>%
  inner_join(i_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(i_IVT_filtered_9_DMSO, by = c("range_start", "range_end"))

A_merged_df_9 <- A_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
i_merged_df_9 <- i_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
A_merged_df_9 <- A_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)
i_merged_df_9 <- i_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

#write.table(A_merged_df_9, "HD10.6_DMSO_STM_A_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(i_merged_df_9, "HD10.6_DMSO_STM_inosine_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
A_plot_long_9 <- melt(A_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")
i_plot_long_9 <- melt(i_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")

unique(A_plot_long_9$Sample)
unique(i_plot_long_9$Sample)

custom_colors_A <- c("black", "black","grey52")
custom_colors_i <- c("purple", "purple", "grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
A_9 <- ggplot(A_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
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
i_9 <-ggplot(i_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_i) +  # Apply custom colors
  labs(x = "modification probability", y = "Inosine fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.02)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
i_9
#m6A from 0.9 to 1.0 probability
i_plot_long_9_filtered <- i_plot_long_9 %>% filter(range_start >= 0.95)

i_9_filtered <-ggplot(i_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_i) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "modification probability", y = "Inosine fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.02)) +  # Set limits and format y-axis
  scale_x_continuous(breaks = seq(0.95, 1, by = 0.01), expand = c(0.003, 0.004)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
i_9_filtered
#combine the plots
combined_plot_9 <- (A_9 | i_9 | i_9_filtered)
combined_plot_9

pdf(file = "D) HD10.6, A&inosine in DMSO&STM in tx&IVT,0.9.0,Profile Plot.pdf", width = 20, height = 5)
combined_plot_9
dev.off()

#DMSO and STM - 0.9.0 - m5C
#Generating the merged dataframe with ranges and fractions
library(dplyr)
C_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "C")
m_tx_filtered_9_DMSO <- HD10_6_DMSO_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("m")& primary_base == "C")
C_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "C")
m_IVT_filtered_9_DMSO <- NHDF_polyA_TSO_IVT.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("m")& primary_base == "C")
C_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("-")& primary_base == "C")
m_tx_filtered_9_STM <- HD10_6_STM2457_48h_1.sup.allMods.trimAdapters.dorado.0.9.0_gencode_v47_probabilities %>% filter (code %in% c("m")& primary_base == "C")


C_merged_df_9 <- C_tx_filtered_9_DMSO %>%
  inner_join(C_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(C_IVT_filtered_9_DMSO, by = c("range_start", "range_end"))

m_merged_df_9 <- m_tx_filtered_9_DMSO %>%
  inner_join(m_tx_filtered_9_STM, by = c("range_start", "range_end")) %>%
  inner_join(m_IVT_filtered_9_DMSO, by = c("range_start", "range_end")) 

C_merged_df_9 <- C_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
m_merged_df_9 <- m_merged_df_9 %>%
  rename(
    DMSO = frac.x,
    STM2457 = frac.y,
    IVT = frac
  )
C_merged_df_9 <- C_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)
m_merged_df_9 <- m_merged_df_9 %>% select(range_start, DMSO, STM2457, IVT)

write.table(C_merged_df_9, "HD10.6_DMSO_STM_C_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)
write.table(m_merged_df_9, "HD10.6_DMSO_STM_m_tx&IVT_0.9.0_merged_df.txt", sep = "\t", row.names = FALSE)

#plotting the probabilities
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Reshape the data to long format for ggplot2
C_plot_long_9 <- melt(C_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")
m_plot_long_9 <- melt(m_merged_df_9, id.vars = "range_start", variable.name = "Sample", value.name = "Value")


custom_colors_C <- c("black", "black","grey52")
custom_colors_m <- c("springgreen4", "springgreen4", "grey52")

custom_linetypes <- c("solid", "dashed", "solid")

# Create the plots
C_9 <- ggplot(C_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_C) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "base probability", y = "C fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
C_9
m_9 <-ggplot(m_plot_long_9, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_linetype_manual(values = custom_linetypes) +
  scale_color_manual(values = custom_colors_m) +  # Apply custom colors
  labs(x = "modification probability", y = "m5C fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.02)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m_9
#m5C from 0.9 to 1.0 probability
m_plot_long_9_filtered <- m_plot_long_9 %>% filter(range_start >= 0.95)

m_9_filtered <-ggplot(m_plot_long_9_filtered, aes_string(x = "range_start", y = "Value", color = "Sample", group = "Sample", linetype = "Sample")) +
  geom_line(size = 1) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors_m) +  # Apply custom colors
  scale_linetype_manual(values = custom_linetypes) +
  labs(x = "modification probability", y = "m5C fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.01)) +  # Set limits and format y-axis
  scale_x_continuous(breaks = seq(0.95, 1, by = 0.01), expand = c(0.003, 0.004)) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 26, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness
m_9_filtered
#combine the plots
combined_plot_9 <- (C_9 | m_9 | m_9_filtered)
combined_plot_9

pdf(file = "E) HD10.6, C&m5C in DMSO&STM in tx&IVT,0.9.0,Profile Plot.pdf", width = 20, height = 5)
combined_plot_9
dev.off()




