library(data.table)
library(grid)
library(extrafont)
library(showtext)
library(sysfonts)
library(grid)
font_import(pattern = "calibri", prompt = FALSE)
loadfonts(device = "pdf")

dorado_NHDF <- fread("dorado_m6A_DRACH_motifs.final.bed")
GLORI_NHDF  <- fread("GLORI_NHDF_motif_DRACH.bed")

dorado_filtered  <- dorado_NHDF[percent_modified_DMSO > percent_modified_STM]
GLORI_NHDF_DRACH <- GLORI_NHDF[(DRACH == TRUE | DRACH == "TRUE") & m6A_level_DMSO >= 0.1]

stopifnot(all(c("chr","start") %in% names(dorado_filtered)),
          all(c("chr","start") %in% names(GLORI_NHDF_DRACH)))

dorado_filtered[, pair := paste(chr, start, sep=":")]
GLORI_NHDF_DRACH[, pair := paste(chr, start, sep=":")]

set_d <- unique(na.omit(dorado_filtered$pair))  
set_g <- unique(na.omit(GLORI_NHDF_DRACH$pair))     

n_d        <- length(set_d)
n_g        <- length(set_g)
n_overlap  <- length(intersect(set_d, set_g))
pct_of_glori <- if (n_g > 0) (100 * n_overlap / n_g) else 0

r_d <- 0.45    # outer Dorado circle (fixed)
r_g <- 0.21    # <<< manually set inner GLORI circle radius (adjust per plot!)

cairo_pdf("Dorado_GLORI_manual_10.pdf", width = 7, height = 8)
grid.newpage()

grid.text(expression("m"^6*"A/A site ratio ≥ 10%"),
          x = 0.5, y = 0.91,              
          just = "top",                   
          gp = gpar(fontsize = 30, family = "Calibri"))

# Outer Dorado
grid.circle(x = 0.5, y = 0.45, r = r_d,
            gp = gpar(fill = "#F59B5B", col = "#F59B5B"))

# Inner GLORI (size = r_g you set above)
grid.circle(x = 0.5, y = 0.45, r = r_g,
            gp = gpar(fill = "#AFEEEE", col = "#AFEEEE"))

# Labels
grid.text("Overlapping\nwith Dorado", x = 0.5, y = 0.54, gp = gpar(fontsize = 27, family = "Calibri"))
grid.text(format(n_overlap, big.mark=","), x = 0.5, y = 0.44, gp = gpar(fontface = "bold", fontsize = 27, family = "Calibri"))
grid.text(sprintf("(%.1f%% of GLORI)", pct_of_glori), x = 0.5, y = 0.38, gp = gpar(fontsize = 23, family = "Calibri"))
grid.text("Dorado DRACH",
          x = 0.5, y = 0.74,
          gp = gpar(fontsize = 24, family = "Calibri"))
          
grid.text(expression("72,372 m"^6*"A sites"),
          x = 0.5, y = 0.68,
          gp = gpar(fontsize = 24, family = "Calibri"))
dev.off()

