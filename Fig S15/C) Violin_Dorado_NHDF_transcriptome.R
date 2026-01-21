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

allContext <- read_delim("NHDF_transcriptome_m6A_all_context_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  # Adjust delimiter if needed
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
  labs(title = "Dorado all-context\n transcriptome NHDF",
       x = expression(paste("91,279 ", m^{6}, "A sites")),
       y = expression("m"^{6}*"A/A site ratio (%)") 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_allContext$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
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
ggsave("violin_trx_NHDFs_m6A_all_context_high-confidence_sites.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed


##### Plot DRACH model

DRACH_model <- read_delim("NHDF_transscriptome_m6A_DRACH_filtered_higher_percent_modified_DMSO.txt", delim = "\t")  # Adjust delimiter if needed
sorted_DRACH_model <- DRACH_model[order(-DRACH_model[[2]]), ]

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
  labs(title = "Dorado DRACH\n transcriptome",
       x = expression(paste("82,325 ", m^{6}, "A sites")),
       y = expression("m"^{6}*"A/A site ratio (%)") 
  ) +
  scale_y_continuous(breaks = seq(0, max(long_sorted_DRACH_model$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
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

ggsave("violin_trx_NHDFs_m6A_DRACH_high-confidence_sites.pdf", plot = p, width = 4, height = 6) 

