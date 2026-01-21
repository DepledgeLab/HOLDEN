library(ggplot2)
library(VennDiagram)
library(readr)
library(data.table)
library(grid)
library(dplyr)
library(ggpointdensity)
library(RColorBrewer)

NHDF_DRACH <- fread("hg38_NHDF_DRACH_highConfidence_with_gene.txt")
colnames(NHDF_DRACH)

HD106_DRACH <- fread("hg38_HD106_DRACH_highConfidence_with_gene.txt")
colnames(HD106_DRACH)
HD106_DRACH$chrom <- sub("^chr", "", HD106_DRACH$chrom)

colnames(NHDF_DRACH)[colnames(NHDF_DRACH) == "chr"]   <- "chrom"
colnames(NHDF_DRACH)[colnames(NHDF_DRACH) == "start"] <- "start_position"
colnames(NHDF_DRACH)[colnames(NHDF_DRACH) == "end"]   <- "end_position"

working_file1 <- merge(NHDF_DRACH, HD106_DRACH, by = c("chrom", "start_position", "end_position"), all=TRUE, sep = "\t")

keys <- c("chrom", "start_position", "end_position")
overlap <- inner_join(NHDF_DRACH, HD106_DRACH, by = keys, suffix = c(".NHDF", ".HD106"))
nhdf_only <- anti_join(NHDF_DRACH, HD106_DRACH, by = keys)
hd106_only <- anti_join(HD106_DRACH, NHDF_DRACH, by = keys)

c(
  NHDF_only = nrow(nhdf_only),
  HD106_only = nrow(hd106_only),
  overlapping = nrow(overlap)
)

colnames(working_file1) <- gsub("(.*)\\.x$", "\\1_NHDF", colnames(working_file1))
colnames(working_file1) <- gsub("(.*)\\.y$", "\\1_HD106", colnames(working_file1))

overlap <- inner_join(NHDF_DRACH, HD106_DRACH, by = keys,
                      suffix = c("-NHDF", "-HD106"))

write.table(nhdf_only,  "NHDF_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(hd106_only, "HD106_only.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(overlap,    "NHDF_HD106_overlap.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


data <- read.table("NHDF_HD106_overlap.txt",
                   header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

data <- na.omit(data)
data <- data[data$gene.NHDF != "" & data$gene.HD106 != "", ]

correlation <- cor(data$`percent_modified_DMSO.NHDF`,
                   data$`percent_modified_DMSO.HD106`,
                   use = "complete.obs")

n_points <- sum(complete.cases(data$`percent_modified_DMSO.NHDF`,
                               data$`percent_modified_DMSO.HD106`))

library(showtext); library(sysfonts)
if (.Platform$OS.type == "windows") {
  font_add("Calibri", regular = "C:/Windows/Fonts/calibri.ttf")
} else if (Sys.info()[["sysname"]] == "Darwin") {
  if (file.exists("/Library/Fonts/Calibri.ttf")) {
    font_add("Calibri", regular = "/Library/Fonts/Calibri.ttf")
  } else if (file.exists("/Library/Fonts/Microsoft/Calibri.ttf")) {
    font_add("Calibri", regular = "/Library/Fonts/Microsoft/Calibri.ttf")
  }
} else {
  font_add("Calibri", regular = "~/.local/share/fonts/calibri.ttf")
}
showtext_auto()

x_name <- if ("percent_modified_DMSO.NHDF" %in% names(data)) "percent_modified_DMSO.NHDF" else "percent_modified_DMSO-NHDF"
y_name <- if ("percent_modified_DMSO.HD106" %in% names(data)) "percent_modified_DMSO.HD106" else "percent_modified_DMSO-HD106"

df <- data[complete.cases(data[[x_name]], data[[y_name]]), ]

r <- cor(df[[x_name]], df[[y_name]], method = "pearson", use = "complete.obs")  

p <- ggplot(df, aes(x = .data[[x_name]], y = .data[[y_name]])) +
  geom_pointdensity(size = 0.6, adjust = 1, method = "neighbors") +
  scale_color_gradientn(
    colours = c("darkblue", "turquoise", "#fee08b", "#fdae61", "#d73027"),
    values  = scales::rescale(c(0, 0.35, 0.6, 0.8, 1)),
    name    = "n_neighbors"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(10, 100), ylim = c(10, 100), expand = FALSE) +
  scale_x_continuous(breaks = seq(10, 100, by = 10)) +
  scale_y_continuous(breaks = seq(10, 100, by = 10)) +
  labs(
    title = expression("38,217 shared m"^6*"A sites"),
    x = expression("NHDF " * "m"^6 * "A level (%)"),
    y = expression("HD10.6 " * "m"^6 * "A level (%)"),
    subtitle = paste0("Pearson r = ", round(r, 4))  
  ) +
  theme_minimal(base_family = "Calibri", base_size = 20) + 
  theme(
    plot.title      = element_text(face = "bold", size = 40),  
    plot.subtitle   = element_text(size = 30),                 
    axis.title.x    = element_text(size = 35, face = "bold"),   
    axis.title.y    = element_text(size = 35, face = "bold"),   
    axis.text.x     = element_text(size = 30),                  
    axis.text.y     = element_text(size = 30),
    legend.title    = element_text(size = 30, face = "bold"),
    legend.text     = element_text(size = 30),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.4, "inch"),
    aspect.ratio = 1.3)

ggsave("NHDF_HD106_neighbors_style.pdf", p, width = 12, height = 10, device = cairo_pdf)
print(p)

