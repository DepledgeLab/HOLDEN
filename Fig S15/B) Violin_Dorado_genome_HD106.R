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

genome_all_context <- read_delim("violin_HD106_m6A_all_context_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  

long_genome_all_context <- genome_all_context %>% 
  pivot_longer(cols = c("DMSO", "STM"), 
               names_to = "Sample",
               values_to = "Value")
str(long_genome_all_context)
head(long_genome_all_context)

long_genome_all_context <- long_genome_all_context %>%
  mutate(Sample = recode(Sample, 
                         "DMSO" = "DMSO", 
                         "STM" = "STM2457"))

p <- ggplot(long_genome_all_context, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15) +  
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  
  scale_fill_manual(values = c("#2874A6", "#CFE2F3")) + 
  labs(title = "Dorado all-context\n (all motifs)",
       x = expression(paste("77,728 ", m^{6}, "A sites")), family = "Calibri",
       y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri" 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_genome_all_context$Value, na.rm = TRUE), by = 10)) +  
  theme_minimal(base_size = 17) + 
  theme(panel.background = element_rect(fill = "white"),  
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none",
        plot.title = element_text(
          size = 24,
          family = "Calibri",
          face = "bold",
          hjust = 0.5   # <-- centers the title
        )
  )

print(p)

ggsave("violin_hg38_HD106_m6A_all_context_high-confidence_sites.pdf", plot = p, width = 4, height = 6)  

##### Plot only DRACH motifs

allContext_DRACH <- read_delim("HD106_m6A_allcontext_filtered_motif_final.bed", delim = "\t", col_names = TRUE)
drach_pattern <- "^[AGU][AG]AC[ACU]$"

drach_matches <- allContext_DRACH %>%
  filter(grepl(drach_pattern, RNA_sequence))

drach_count <- nrow(drach_matches)
cat("Number of DRACH motif matches:", drach_count, "\n")

write_delim(drach_matches, "HD106_allContext_DRACH_motifs_extracted.bed", delim = "\t")
allContext_DRACH <- read_delim("HD106_allContext_DRACH_motifs_extracted.bed", delim = "\t", col_names = TRUE)

long_allContext_DRACH <- allContext_DRACH %>%
  pivot_longer(cols = c("percent_modified_DMSO", "percent_modified_STM"),  
               names_to = "Sample",
               values_to = "Value")
str(long_allContext_DRACH)
head(long_allContext_DRACH)

long_allContext_DRACH <- long_allContext_DRACH %>%
  mutate(Sample = recode(Sample, 
                         "percent_modified_DMSO" = "DMSO", 
                         "percent_modified_STM" = "STM2457"))

p <- ggplot(long_allContext_DRACH, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15) +  
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  
  scale_fill_manual(values = c("#2874A6", "#CFE2F3")) + 
  labs(title = "Dorado all-context\n (DRACH only)",
       x = expression(paste("67,068 ", m^{6}, "A sites")), family = "Calibri",
       y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri" 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_allContext_DRACH$Value, na.rm = TRUE), by = 10)) + 
  theme_minimal(base_size = 17) +  
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none",
        plot.title = element_text(
          size = 24,
          family = "Calibri",
          face = "bold",
          hjust = 0.5   # <-- centers the title
        )
  )

print(p)
ggsave("violin_hg38_HD106_m6A_all_context_high-confidence_sites_DRACH_only.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed




##### Plot DRACH model

DRACH_model <- read_delim("violin_HD106_m6A_DRACH_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  
sorted_DRACH_model <- DRACH_model[order(-DRACH_model[[2]]), ]

long_sorted_DRACH_model <- sorted_DRACH_model %>%
  pivot_longer(cols = c("percent_modified_DMSO", "percent_modified_STM"), 
               names_to = "Sample",
               values_to = "Value")
str(long_sorted_DRACH_model)
head(long_sorted_DRACH_model)

long_sorted_DRACH_model <- long_sorted_DRACH_model %>%
  mutate(Sample = recode(Sample, 
                         "percent_modified_DMSO" = "DMSO", 
                         "percent_modified_STM" = "STM2457"))

p <- ggplot(long_sorted_DRACH_model, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) + 
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  
  scale_fill_manual(values = c("#2874A6", "#CFE2F3")) + 
  labs(
    title = "Dorado DRACH\nmodel",
    x = expression(paste("67,222 ", m^{6}, "A sites")),
    y = expression("m"^{6}*"A/A site ratio (%)")
  ) +
  scale_y_continuous(breaks = seq(0, max(long_sorted_DRACH_model$Value, na.rm = TRUE), by = 10)) +  
  theme_minimal(base_size = 17) +  
  theme(
    panel.background = element_rect(fill = "white"),  
    panel.grid.major.y = element_line(color = "grey", size = 0.05),
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(color = "black", size = 20, family = "Calibri"),
    axis.text.y = element_text(color = "black", size = 20, family = "Calibri"),
    axis.title.x = element_text(size = 24, family = "Calibri"), 
    axis.title.y = element_text(size = 24, family = "Calibri"),
    legend.position = "none",
    plot.title = element_text(      # ← keep only this one
      size = 24,
      family = "Calibri",
      face = "bold",
      hjust = 0.5
    )
  )

print(p)
ggsave("violin_hg38_HD106_m6A_DRACH_high-confidence_sites.pdf", plot = p, width = 4, height = 6)  

