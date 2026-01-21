library(UpSetR)
library(ggplot2)  
library(scales)   
library(data.table)

dorado <- fread(
  "NHDF_genome_DRACH_motif_annotated.txt")

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

library(UpSetR)
library(showtext)
library(grid)
library(data.table)

set_order <- c("GLORI", "GLORI_DRACH", "Dorado_DRACH")

all_sites <- unique(c(set_glori, set_drach, set_dorado))
upset_dt <- data.table(
  site_id = all_sites,
  GLORI        = as.integer(all_sites %in% set_glori),
  GLORI_DRACH  = as.integer(all_sites %in% set_drach),
  Dorado_DRACH = as.integer(all_sites %in% set_dorado))

font_add(family = "Calibri",
         regular    = "C:/Windows/Fonts/calibri.ttf",
         bold       = "C:/Windows/Fonts/calibrib.ttf",
         italic     = "C:/Windows/Fonts/calibrii.ttf",
         bolditalic = "C:/Windows/Fonts/calibriz.ttf")
showtext_auto()
pdf.options(useDingbats = FALSE)

pdf("UpSet_GLORI_Dorado_final_singlepage.pdf", width = 35, height = 12, onefile = FALSE)
showtext_begin()

options(scipen = 999)
par(family = "Calibri", cex = 5, mar = c(8, 8, 8, 8))

upset(
  upset_dt,
  sets = set_order,
  keep.order = TRUE,
  sets.bar.color = "grey60",
  main.bar.color = "grey30",
  point.size = 5,
  line.size = 2,
  order.by = "freq",
  empty.intersections = "on",
  text.scale = c(5, 5, 5, 5, 5, 5),
  mainbar.y.label = "Intersection size",
  sets.x.label = "Set size",
  mb.ratio =  c(0.65, 0.35))

# --- Robust relabeling (handles all UpSetR versions) ---
grid.force()

rename_map <- list(
  "Dorado_DRACH" = "Dorado",
  "GLORI_DRACH"  = "GLORI (DRACH)",
  "GLORI"        = "GLORI"
)

# Find all grobs (axes and texts)
all_grobs <- grid.ls(print = FALSE)$name
candidate_grobs <- grep("(axis|text)", all_grobs, value = TRUE)

# Search for and replace the target labels
for (g in candidate_grobs) {
  grob_obj <- try(grid.get(g), silent = TRUE)
  if (inherits(grob_obj, "grob") && !is.null(grob_obj$label)) {
    labels_vec <- as.character(grob_obj$label)
    new_labels_vec <- vapply(labels_vec, function(lbl) {
      if (lbl %in% names(rename_map)) rename_map[[lbl]] else lbl
    }, character(1))
    if (!identical(labels_vec, new_labels_vec)) {
      grid.edit(g, label = new_labels_vec,
                gp = gpar(fontfamily = "Calibri", fontsize = 40, col = "black"))
    }
  }
}

grid.refresh()
showtext_end()
dev.off()

