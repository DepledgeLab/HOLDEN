library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(showtext)

font_add("Calibri", regular = "C:/Windows/Fonts/calibri.ttf") 
showtext_auto()

data <- read_delim("NHDF_m6A_DRACH_motifs.final.bed", delim = "\t", col_names = TRUE)

drach_pattern <- "^[AGU][AG]AC[ACU]$"
data <- data %>%
  mutate(DRACH_status = ifelse(grepl(drach_pattern, RNA_sequence), RNA_sequence, "non-DRACH"))

motif_counts <- data %>%
  count(DRACH_status) %>%
  arrange(desc(n))

motif_levels <- c(setdiff(motif_counts$DRACH_status, "non-DRACH"), "non-DRACH")
data$DRACH_status <- factor(data$DRACH_status, levels = motif_levels)
motif_counts$DRACH_status <- factor(motif_counts$DRACH_status, levels = motif_levels)

# Create color palette
n_drach <- length(motif_levels) - 1
color_vector <- c(colorRampPalette(c("#EB984E", "#F8C471", "#FAD7A0"))(n_drach), "gray60")
names(color_vector) <- motif_levels

pdf("NHDF_m6A_violin_plots_motifs_DRACH_model.pdf", width = 8, height = 5, useDingbats = FALSE)

ggplot(data, aes(x = DRACH_status, y = percent_modified_DMSO, fill = DRACH_status)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, fill = "white", color = "black", alpha = 0.1) +
  scale_fill_manual(values = color_vector) +
  theme_minimal(base_size = 14, base_family = "Calibri") +
  labs(
    x = "5-mer motif",
    y = expression("m"^6*"A/A site ratio (%)"),
    #title = expression("NHDF - Dorado m"^6*"A DRACH model")
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20, family = "Calibri"),
    axis.text.y = element_text(size = 20, family = "Calibri"),
    axis.title.x = element_text(size = 22, family = "Calibri"),
    axis.title.y = element_text(size = 22, family = "Calibri"),
    plot.title = element_text(size = 24, hjust = 0.5, family = "Calibri"),
    legend.position = "none")

dev.off()


#Barplot 
major_breaks <- 10^seq(0, 5)                
minor_breaks <- as.numeric(outer(1:9, 10^(0:5)))  

pdf("NHDF_m6A_motif_counts_barplot.pdf", width = 8, height = 6, useDingbats = FALSE)

ggplot(motif_counts, aes(x = DRACH_status, y = n, fill = DRACH_status)) +
  geom_bar(stat = "identity", width = 0.7, color = NA) +
  scale_fill_manual(values = color_vector) +
  scale_y_log10(
    breaks = major_breaks,
    minor_breaks = minor_breaks,
    labels = comma
  ) +
  theme_minimal(base_size = 14, base_family = "Calibri") +
  labs(
    x = "5-mer motif",
    y = expression("Number of m"^6*"A sites"),
    title = expression("NHDF - Dorado m"^6*"A DRACH model")
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22, family = "Calibri"),
    axis.text.y = element_text(size = 22, family = "Calibri"),
    axis.title.x = element_text(size = 24, family = "Calibri"),
    axis.title.y = element_text(size = 24, family = "Calibri"),
    plot.title = element_text(size = 25, hjust = 0.5, family = "Calibri"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray70", size = 0.3),  
    panel.grid.minor.y = element_line(color = "gray90", size = 0.2),   
    panel.grid.major.x = element_line(color = "gray85", size = 0.3),
    panel.grid.minor.x = element_blank())

dev.off()

