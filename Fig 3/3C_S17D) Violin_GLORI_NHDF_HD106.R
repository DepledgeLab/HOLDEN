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

GLORI_all_NHDFs <- read_delim("GLORI_NHDF_motif_DRACH.bed", delim = "\t", show_col_types = FALSE)  

long_GLORI_all_NHDFs <- GLORI_all_NHDFs %>% 
  pivot_longer(cols = c("m6A_level_DMSO", "m6A_level_STM"),  
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Value = Value * 100)

str(long_GLORI_all_NHDFs)
head(long_GLORI_all_NHDFs)

long_GLORI_all_NHDFs <- long_GLORI_all_NHDFs %>%
  mutate(Sample = recode(Sample, 
                         "m6A_level_DMSO" = "DMSO", 
                         "m6A_level_STM" = "STM2457"))

p <- ggplot(long_GLORI_all_NHDFs, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("turquoise", "#AFEEEE")) + 
  labs(
    title = "GLORI\n (all motifs)",
    x = expression("140,"*934~m^{6}*"A sites"), family = "Calibri",
    y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri"      
  ) +
  scale_y_continuous(breaks = seq(0, max(long_GLORI_all_NHDFs$Value, na.rm = TRUE), by = 10)) + 
  theme_minimal(base_size = 17) +  # Base font size
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "grey", size = 0.05), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none" )

print(p)
ggsave("violin_GLORI_NHDF_high-confidence_sites.pdf", plot = p, width =4, height = 6)  



### now only DRACH motifs
GLORI_DRACH_NHDFs <- GLORI_all_NHDFs %>% 
  filter(DRACH == "TRUE")

long_GLORI_DRACH_NHDFs <- GLORI_DRACH_NHDFs %>% 
  pivot_longer(cols = c("m6A_level_DMSO", "m6A_level_STM"),
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Value = Value * 100) %>%
  mutate(Sample = recode(Sample, 
                         "m6A_level_DMSO" = "DMSO", 
                         "m6A_level_STM" = "STM2457"))

p <- ggplot(long_GLORI_DRACH_NHDFs, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) + 
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) + 
  scale_fill_manual(values = c("turquoise", "#AFEEEE")) + 
  labs(
    title = "GLORI\n (DRACH motifs)",
    x = expression("27,"*926~m^{6}*"A sites"), family = "Calibri",      
    y = expression("m"^{6}*"A/A site ratio (%)"), family = "Calibri"    
  ) +
  scale_y_continuous(breaks = seq(0, max(long_GLORI_DRACH_NHDFs$Value, na.rm = TRUE), by = 10)) +  
  theme_minimal(base_size = 17) +  # Base font size
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "grey", size = 0.05),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        legend.position = "none" )

print(p)
ggsave("violin_GLORI_NHDF_DRACH_only.pdf", plot = p, width = 4, height = 6)  


#### Same for HD106

GLORI_all_HD106 <- read_delim("GLORI_HD106_motif_DRACH.bed", delim = "\t", show_col_types = FALSE)  # Adjust delimiter if needed

long_GLORI_all_HD106 <- GLORI_all_HD106 %>% 
  pivot_longer(cols = c("m6A_level_DMSO", "m6A_level_STM"),  # Replace with actual column names
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Value = Value * 100)

str(long_GLORI_all_HD106)
head(long_GLORI_all_HD106)

# Customize sample labels
long_GLORI_all_HD106 <- long_GLORI_all_HD106 %>%
  mutate(Sample = recode(Sample, 
                         "m6A_level_DMSO" = "DMSO", 
                         "m6A_level_STM" = "STM2457"))
# Create the violin plot
p <- ggplot(long_GLORI_all_HD106, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("turquoise", "#AFEEEE")) + 
  labs(title = "HD10.6\nGLORI (all motifs)",
       x = expression("174,"*282~m^{6}*"A sites"),   # or "174,"*282 if that’s the correct count here
       y = expression("m"^{6}*"A/A site ratio (%)")
  ) +
  scale_y_continuous(breaks = seq(0, max(long_GLORI_all_HD106$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
  theme_minimal(base_size = 17) +  # Base font size
  theme(panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  # Major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

print(p)
ggsave("violin_GLORI_HD106_high-confidence_sites.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed


### now only DRACH motifs
GLORI_DRACH_HD10 <- GLORI_all_HD106 %>% 
  filter(DRACH == "TRUE")

long_GLORI_DRACH_HD10 <- GLORI_DRACH_HD10 %>% 
  pivot_longer(cols = c("m6A_level_DMSO", "m6A_level_STM"),
               names_to = "Sample",
               values_to = "Value") %>%
  mutate(Value = Value * 100) %>%
  mutate(Sample = recode(Sample, 
                         "m6A_level_DMSO" = "DMSO", 
                         "m6A_level_STM" = "STM2457"))

p <- ggplot(long_GLORI_DRACH_HD10, aes(x = Sample, y = Value, fill = Sample)) +
  geom_violin(trim = TRUE, scale = "width", width = 1, alpha = 1) +  # Adjust width here
  geom_boxplot(width = 0.05, color = "black", alpha = 1, size = 0.15, outlier.size = 0.1) +  # Thinner boxplot lines
  geom_vline(xintercept = 1.5, color = "black", size = 0.05) +  # Vertical line between samples
  scale_fill_manual(values = c("turquoise", "#AFEEEE")) + 
  labs(title = "HD10.6\nGLORI (DRACH motifs)",
       x = expression("27,"*965~m^{6}*"A sites"),   
       y = expression("m"^{6}*"A/A site ratio (%)")
  ) +
  scale_y_continuous(breaks = seq(0, max(long_GLORI_DRACH_HD10$Value, na.rm = TRUE), by = 10)) +  # Set y-axis breaks to intervals of 10
  theme_minimal(base_size = 17) +  # Base font size
  theme(panel.background = element_rect(fill = "white"),  # Set background to white
        panel.grid.major.y = element_line(color = "grey", size = 0.05),  # Major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        axis.title.x = element_text(size = 24, family = "Calibri"), 
        axis.title.y = element_text(size = 24, family = "Calibri"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

print(p)
ggsave("violin_GLORI_HD106_DRACH_only.pdf", plot = p, width = 4, height = 6)  # Adjust dimensions as needed



