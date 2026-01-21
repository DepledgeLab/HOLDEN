library(UpSetR)
library(ggplot2)  
library(scales)   
library(data.table)
library(showtext)
library(grid)

dorado <- fread ("NHDF_genome_DRACH_motif_annotated.txt")

glori <- fread(
  "GLORI_NHDF_motif_DRACH.bed",
  col.names = c("chr", "start", "end", "Gene", "ENSG", "strand",
                "coverage_DMSO", "m6A_level_DMSO", "P_value_DMSO",
                "coverage_STM", "m6A_level_STM", "P_value_STM",
                "RNA_sequence", "DRACH"))

dorado_filtered <- dorado[percent_modified_DMSO >= percent_modified_STM]
glori_DRACH <- glori[DRACH == TRUE | DRACH == "TRUE" | DRACH == "True"]
glori_filtered <- glori[!is.na(m6A_level_DMSO) & !is.na(m6A_level_STM) &
                          m6A_level_DMSO >= m6A_level_STM]

glori_filtered[, site_id := paste(chr, start, end, sep = "_")]
glori_DRACH[,    site_id := paste(chr, start, end, sep = "_")]
dorado_filtered[, site_id := paste(chrom, start, end, sep = "_")]

set_glori   <- unique(glori_filtered$site_id)
set_drach   <- unique(glori_DRACH$site_id)
set_dorado  <- unique(dorado_filtered$site_id)

n_glori   <- length(set_glori)
n_drach   <- length(set_drach)
n_dorado  <- length(set_dorado)

n_glori_drach  <- length(intersect(set_glori, set_drach))
n_glori_dorado <- length(intersect(set_glori, set_dorado))
n_drach_dorado <- length(intersect(set_drach, set_dorado))

n_all_three <- length(Reduce(intersect, list(set_glori, set_drach, set_dorado)))

cat("=== m⁶A site counts ===\n")
cat("GLORI total:             ", n_glori, "\n")
cat("GLORI (DRACH) total:     ", n_drach, "\n")
cat("Dorado (DRACH model):    ", n_dorado, "\n\n")

cat("GLORI ∩ GLORI(DRACH):    ", n_glori_drach, "\n")
cat("GLORI ∩ Dorado:          ", n_glori_dorado, "\n")
cat("GLORI(DRACH) ∩ Dorado:   ", n_drach_dorado, "\n")
cat("All three:               ", n_all_three, "\n")



font_add(family = "Calibri",
         regular    = "C:/Windows/Fonts/calibri.ttf",
         bold       = "C:/Windows/Fonts/calibrib.ttf",
         italic     = "C:/Windows/Fonts/calibrii.ttf",
         bolditalic = "C:/Windows/Fonts/calibriz.ttf")
showtext_auto()

m6A_overlap <- c(
  "GLORI (all motifs)" = 140934,
  "GLORI (DRACH motifs)" = 27926,
  "Dorado (DRACH model)" = 72372,
  "GLORI (all motifs)&GLORI (DRACH motifs)" = 27926,
  "GLORI (all motifs)&Dorado (DRACH model)" = 13631,
  "GLORI (DRACH motifs)&Dorado (DRACH model)" = 13614,
  "GLORI (all motifs)&GLORI (DRACH motifs)&Dorado (DRACH model)" = 13614)

input_expr <- UpSetR::fromExpression(m6A_overlap)

pdf("UpSet_GLORI_DORADO.pdf", width = 22, height = 8)
showtext_begin()

UpSetR::upset(
  input_expr,
  order.by = "freq",
  decreasing = TRUE,
  mb.ratio = c(0.65, 0.35),
  text.scale = 4.2,         # large text across the figure
  point.size = 5,           # larger intersection dots
  line.size = 2,            # thicker connecting lines
  show.numbers = FALSE,
  sets.x.label = "",
  sets.bar.color = "gray27",
  main.bar.color = "darkgrey",
  matrix.color = "gray27")

grid.force()
grid.text(
  label = prettyNum(m6A_overlap, big.mark = ",", preserve.width = "none"),
  x = unit(seq_along(m6A_overlap), "native"),
  y = unit(m6A_overlap + max(m6A_overlap) * 0.04, "native"),
  gp = gpar(fontsize = 50, fontfamily = "Calibri", fontface = "bold", col = "black"))

showtext_end()
dev.off()

