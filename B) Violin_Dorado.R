library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(data.table)
library(showtext)

font_add("Calibri", regular = "C:/Windows/Fonts/calibri.ttf")
showtext_auto()

allContext <- read_delim("violin_hg38_NHDFs_m6A_all_context_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  # Adjust delimiter if needed
sorted_allContext <- allContext[order(-allContext[[2]]), ]
print(sorted_allContext)

long_allContext <- sorted_allContext %>%
  pivot_longer(cols = c("percent_modified_DMSO", "percent_modified_STM"),  # Replace with actual column names
               names_to = "Sample",
               values_to = "Value")
str(long_allContext)
head(long_allContext)

long_allContext <- long_allContext %>%
  mutate(Sample = recode(Sample, 
                         "percent_modified_DMSO" = "DMSO", 
                         "percent_modified_STM" = "STM2457"))

p <- ggplot(long_allContext, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("#E56A00", "#F5B59B")) + 
  labs(title = "Dorado all-context\n (all motifs)",
       x = expression(paste("81,012 ", m^{6}, "A sites")), family = "Calibri",
       y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri" 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_allContext$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
  theme_minimal(base_size = 17) +  # Base font size
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  # Major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"), 
        legend.position = "none" )

print(p)
ggsave("violin_hg38_NHDFs_m6A_all_context_high-confidence_sites.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed


##### Plot only DRACH motifs
allContext_DRACH <- read_delim("m6A_allcontext_motifs.final.bed", delim = "\t", col_names = TRUE)
drach_pattern <- "^[AGU][AG]AC[ACU]$"

drach_matches <- allContext_DRACH %>%
  filter(grepl(drach_pattern, RNA_sequence))

drach_count <- nrow(drach_matches)
cat("Number of DRACH motif matches:", drach_count, "\n")

write_delim(drach_matches, "allContext_DRACH_motifs_extracted.bed", delim = "\t")
allContext_DRACH <- read_delim("allContext_DRACH_motifs_extracted.bed", delim = "\t", col_names = TRUE)

long_allContext_DRACH <- allContext_DRACH %>%
  pivot_longer(cols = c("percent_modified_DMSO", "percent_modified_STM"),  # Replace with actual column names
               names_to = "Sample",
               values_to = "Value")
str(long_allContext_DRACH)
head(long_allContext_DRACH)

long_allContext_DRACH <- long_allContext_DRACH %>%
  mutate(Sample = recode(Sample, 
                         "percent_modified_DMSO" = "DMSO", 
                         "percent_modified_STM" = "STM2457"))


p <- ggplot(long_allContext_DRACH, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("#E56A00", "#F5B59B")) + 
  labs(title = "Dorado all-context\n (DRACH only)",
       x = expression(paste("72,047 ", m^{6}, "A sites")), family = "Calibri",
       y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri" 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_allContext_DRACH$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
  theme_minimal(base_size = 17) +  # Base font size
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  # Major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none" )

print(p)

ggsave("violin_hg38_NHDFs_m6A_allContext_high-confidence_sites_onlyDRACHmotifs.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed




##### Plot DRACH model

DRACH_model <- read_delim("violin_hg38_NHDFs_m6A_DRACH_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  # Adjust delimiter if needed
sorted_DRACH_model <- DRACH_model[order(-data[[2]]), ]

# Reshape the data to long format
long_sorted_DRACH_model <- sorted_DRACH_model %>%
  pivot_longer(cols = c("percent_modified_DMSO", "percent_modified_STM"),  # Replace with actual column names
               names_to = "Sample",
               values_to = "Value")
str(long_sorted_DRACH_model)
head(long_sorted_DRACH_model)

# Calculate the average for each sample
averages <- long_sorted_DRACH_model %>%
  group_by(Sample) %>%
  summarise(Average = mean(Value, na.rm = TRUE))

# Customize sample labels
long_sorted_DRACH_model <- long_sorted_DRACH_model %>%
  mutate(Sample = recode(Sample, 
                         "percent_modified_DMSO" = "DMSO", 
                         "percent_modified_STM" = "STM2457"))

p <- ggplot(long_sorted_DRACH_model, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("#E56A00", "#F5B59B")) + 
  labs(title = "Dorado DRACH\n model",
       x = expression(paste("72,372 ", m^{6}, "A sites")), family = "Calibri",
       y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri" 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_sorted_DRACH_model$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
  theme_minimal(base_size = 17) +  # Base font size
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  # Major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none" )

print(p)

ggsave("violin_hg38_NHDFs_m6A_DRACH_high-confidence_sites.pdf", plot = p, width = 4, height = 6) 

