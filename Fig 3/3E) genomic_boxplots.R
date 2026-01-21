library(Gviz)
library(GenomicFeatures)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(grid)
library(showtext)
library(scales)
library(ggplot2)
library(scales)
library(ggnewscale)

font_add(family = "Calibri",
         regular = "C:/Windows/Fonts/calibri.ttf",
         bold    = "C:/Windows/Fonts/calibrib.ttf",
         italic  = "C:/Windows/Fonts/calibrii.ttf",
         bolditalic = "C:/Windows/Fonts/calibriz.ttf")
showtext_auto()

gff         <- "SPEN.gff3"
dorado_path <- "dorado_genome_csv.txt"   
glori_path  <- "glori_genome_csv.txt"   
target_tx   <- "ENST00000375759.8"

box_bp     <- 30L
gray_line  <- "#BFBFBF"
dorado_pal <- c("#FFF2E0", "#FFDAB3", "#E56A00", "#7A3300")
glori_pal  <- c("#5CB6D3", "#2C89A7", "#0E5A71")
ST_MIN <- 10; ST_MAX <- 100

# ---- Build TxDb and extract ONLY the target isoform ----
txdb <- makeTxDbFromGFF(gff, format = "gff3")
try({ GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC" }, silent = TRUE)

ex_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
nm  <- names(ex_by_tx)
hit <- nm == target_tx
if (!any(hit)) hit <- sub("\\.\\d+$", "", nm) == sub("\\.\\d+$", "", target_tx)
stopifnot(any(hit))
ex_sel <- ex_by_tx[which(hit)[1]]
names(ex_sel) <- target_tx

rng  <- range(unlist(ex_sel, use.names = FALSE))
chr  <- as.character(seqnames(rng))[1]
from <- max(1L, start(rng) - 300L)
to   <- end(rng) + 300L

# ---- Tracks ----
axisTrack <- GenomeAxisTrack(
  littleTicks = TRUE,
  labelPos = "alternating")

displayPars(axisTrack) <- list(
  fontfamily = "Calibri",
  fontsize = 18,        
  cex = 1.4,            
  lwd = 1.2,            
  col = "#555555",      
  col.axis = "black",  
  col.title = "black")

geneTrack <- GeneRegionTrack(
  ex_sel,
  chromosome = chr, genome = "hg38",
  name = "",
  exonOnly = TRUE,
  collapseTranscripts = TRUE, stacking = "dense",
  fill = "#EAEAEA", col = "#888888")

displayPars(geneTrack) <- list(
  lwd              = 0.6,
  col.line         = gray_line,
  col.arrow        = gray_line,
  background.title = "white",
  col.title        = "black",
  fontfamily           = "Calibri",
  thinBoxFeature   = TRUE,
  arrowHeadWidth   = 6,
  fontsize.group   = 20,
  min.height = 2)

# ---- Helper for box tracks ----
harmonize_chr <- function(x_chr, locus_chr) {
  x_chr     <- as.character(x_chr)
  locus_chr <- as.character(locus_chr)[1]
  locus_has <- grepl("^chr", locus_chr)
  x_has     <- length(x_chr) > 0 && grepl("^chr", x_chr[1])
  if (locus_has && !x_has) paste0("chr", x_chr)
  else if (!locus_has && x_has) sub("^chr", "", x_chr)
  else x_chr}

make_box_track <- function(path, label, palette_cols, locus_chr, locus_from, locus_to) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (!all(c("chr","start","strand","stoich") %in% names(dt))) return(NULL)
  
  dt[, chr    := harmonize_chr(chr, locus_chr)]
  dt[, start  := as.integer(start)]
  dt[, strand := as.character(strand)]
  dt[, stoich := suppressWarnings(as.numeric(stoich))]
  dt <- dt[!is.na(chr) & !is.na(start) & !is.na(stoich)]
  
  if (max(dt$stoich, na.rm = TRUE) <= 1.5) dt[, stoich := stoich * 100]
  dt[, stoich := pmax(ST_MIN, pmin(ST_MAX, stoich))]
  dt <- dt[chr == locus_chr & start >= locus_from & start <= locus_to]
  if (!nrow(dt)) return(NULL)
  
  st <- pmax(locus_from, dt$start - floor(box_bp/2))
  en <- pmin(locus_to,   dt$start + ceiling(box_bp/2))
  s  <- dt$strand; s[!s %in% c("+","-","*")] <- "*"
  
  gr <- GRanges(seqnames = dt$chr, ranges = IRanges(st, en), strand = s)
  gr$stoich <- dt$stoich
  
  pal  <- scales::col_numeric(palette = palette_cols, domain = c(ST_MIN, ST_MAX))
  cols <- pal(gr$stoich)
  
  tr <- AnnotationTrack(
    gr,
    chromosome = locus_chr, genome = "hg38",
    name = label,
    shape = "box",
    fill = cols, col  = cols, lwd = 0.2,
    stacking = "dense",          # ← one-row, no vertical staggering
    thinBoxFeature = TRUE,       # ← optional: slimmer boxes
    cex = 1  )
  
  attr(tr, "pal") <- pal   # store color scale for legend
  tr}

doradoTrack <- make_box_track(dorado_path, "Dorado", dorado_pal, chr, from, to)
gloriTrack  <- make_box_track(glori_path,  "GLORI",  glori_pal,  chr, from, to)
if (!is.null(doradoTrack)) displayPars(doradoTrack)$rotation.title <- 0
if (!is.null(gloriTrack))  displayPars(gloriTrack)$rotation.title  <- 0
if (!is.null(doradoTrack)) {
  displayPars(doradoTrack) <- list(
    cex.title = 2.4,            # ← increase this for larger label
    fontfamily.title = "Calibri",
    fontface.title = "plain",
    rotation.title = 0  )}

if (!is.null(gloriTrack)) {
  displayPars(gloriTrack) <- list(
    cex.title = 2.4,
    fontfamily.title = "Calibri",
    fontface.title = "plain",
    rotation.title = 0  )}

# ---- Combine and plot ----
tracks <- list(axisTrack)
sizes  <- c(0.55)            

if (!is.null(doradoTrack)) { tracks <- c(tracks, list(doradoTrack)); sizes <- c(sizes, 0.80) }
if (!is.null(gloriTrack))  { tracks <- c(tracks,  list(gloriTrack));  sizes <- c(sizes, 0.80) }

keep  <- vapply(tracks, function(x) !is.null(x), logical(1))
tracks <- tracks[keep]; sizes <- sizes[keep]
# ---- Make labels horizontal and larger ----
if (!is.null(doradoTrack)) displayPars(doradoTrack) <- list(
  cex.title = 2,          
  fontfamily = "Calibri",
  rotation.title = 0)
if (!is.null(gloriTrack)) displayPars(gloriTrack) <- list(
  cex.title = 2,
  fontfamily = "Calibri",
  rotation.title = 0)

# ---- Colorbar helper (vertical or horizontal) ----
draw_colorbar <- function(pal_fn, title, x, y,
                          orientation = c("vertical","horizontal"),
                          bar_len = 0.36, bar_thick = 0.022, n = 200,
                          tick_cex = 0.9, title_cex = 0.9,
                          # nudges:
                          tick_dx_left = 0, tick_dx_right = 0, tick_dy = 0, title_dy = 0) {
  orientation <- match.arg(orientation)
  vals <- seq(ST_MIN, ST_MAX, length.out = n); colv <- pal_fn(vals)
  
  if (orientation == "horizontal") {
    bar <- matrix(colv, nrow = 26, ncol = n, byrow = TRUE)
    vp <- grid::viewport(x = grid::unit(x,"npc"), y = grid::unit(y,"npc"),
                         width = grid::unit(bar_len,"npc"), height = grid::unit(bar_thick,"npc"),
                         just = c("left","center"), clip = "off")
    grid::pushViewport(vp)
    grid::grid.raster(bar, x = unit(0,"npc"), y = unit(0,"npc"),
                      width = unit(1,"npc"), height = unit(1,"npc"),
                      just = c("left","bottom"), interpolate = FALSE)
    grid::grid.rect(gp = grid::gpar(fill = NA, col = "black", lwd = 0.4))
    
    # ticks (move with tick_dx_* and tick_dy)
    grid::grid.text(sprintf("%g", ST_MIN), x = unit(0 + tick_dx_left,"npc"),
                    y = unit(-0.35 + tick_dy,"npc"), just = "left",
                    gp = grid::gpar(cex = tick_cex, fontfamily = "Calibri"))
    grid::grid.text(sprintf("%g", ST_MAX), x = unit(1 + tick_dx_right,"npc"),
                    y = unit(-0.35 + tick_dy,"npc"), just = "right",
                    gp = grid::gpar(cex = tick_cex, fontfamily = "Calibri"))
    # title
    grid::grid.text(title, x = unit(0.5,"npc"), y = unit(1.45 + title_dy,"npc"),
                    just = "center", gp = gpar(cex = title_cex, fontfamily = "Calibri"))
    grid::upViewport()
    
  } else {
    bar <- matrix(colv, nrow = n, ncol = 26, byrow = FALSE)[n:1, , drop = FALSE]
    vp <- grid::viewport(x = unit(x,"npc"), y = unit(y,"npc"),
                         width = unit(bar_thick,"npc"), height = unit(bar_len,"npc"),
                         just = c("left","center"), clip = "off")
    grid::pushViewport(vp)
    grid::grid.raster(bar, x = unit(0,"npc"), y = unit(0,"npc"),
                      width = unit(1,"npc"), height = unit(1,"npc"),
                      just = c("left","bottom"), interpolate = FALSE)
    grid::grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.4))
    
    grid::grid.text(sprintf("%g", ST_MAX), x = unit(1.15 + tick_dx_right,"npc"),
                    y = unit(1 + tick_dy,"npc"), just = "left",
                    gp = gpar(cex = tick_cex, fontfamily = "Calibri"))
    grid::grid.text(sprintf("%g", ST_MIN), x = unit(1.15 + tick_dx_right,"npc"),
                    y = unit(0 + tick_dy,"npc"), just = "left",
                    gp = gpar(cex = tick_cex, fontfamily = "Calibri"))
    grid::grid.text(title, x = unit(0.5,"npc"), y = unit(1.16 + title_dy,"npc"),
                    just = "center", gp = gpar(cex = title_cex, fontfamily = "Calibri"))
    grid::upViewport()  }}

if (!is.null(doradoTrack)) displayPars(doradoTrack)$fontfamily <- "Calibri"
if (!is.null(gloriTrack))  displayPars(gloriTrack)$fontfamily  <- "Calibri"
displayPars(axisTrack)$fontfamily <- "Calibri"

if (!is.null(doradoTrack)) {
  dp <- displayPars(doradoTrack)
  dp$cex.title       <- 2
  dp$fontfamily      <- "Calibri"
  dp$rotation.title  <- 0
  dp$fontface.title  <- "1"   # ← make title not bold (1 or "plain")
  dp$fontface        <- "plain"   # ← fallback: some Gviz versions use global fontface
  displayPars(doradoTrack) <- dp}

if (!is.null(gloriTrack)) {
  dp <- displayPars(gloriTrack)
  dp$cex.title       <- 2
  dp$fontfamily      <- "Calibri"
  dp$rotation.title  <- 0
  dp$fontface.title  <- "1"   # ← not bold
  dp$fontface        <- "plain"   # ← fallback
  displayPars(gloriTrack) <- dp}

pdf("SPEN_ENST00000375759.8_genomic.pdf", width = 12, height = 4.5, useDingbats = FALSE)
for (tr in tracks) {
  if (inherits(tr, "Track")) {
    dp <- displayPars(tr)
    dp$fontfamily <- "Calibri"
    dp$fontfamily.title <- "Calibri"
    dp$fontface.title <- "plain"
    dp$rotation.title <- 0
    dp$cex.title <- 2
    displayPars(tr) <- dp  }}
displayPars(axisTrack)$fontfamily <- "Calibri"
displayPars(geneTrack)$fontfamily <- "Calibri"

plotTracks(tracks,
           chromosome = chr, from = from, to = to,
           sizes = sizes,
           background.title = "white", col.title = "black",
           extend.left = 100,
           cex = 1, fontfamily = "Calibri",
           cex.title = 2,
           fontface.title = "plain",
           margin = 28)
displayPars(doradoTrack)$size <- 0.5
displayPars(gloriTrack)$size  <- 0.5

grid.text(expression(italic("SPEN") ~" m"^6*"A sites (genomic alignment)"),
          x = 0.5, y = 0.94,
          gp = gpar(fontsize = 24, fontfamily = "Calibri"))

draw_colorbar(attr(doradoTrack,"pal"),
              expression("Dorado " * m^6 * "A/A (%)"),
              x = 0.5, y = 0.5,
              orientation = "horizontal",
              bar_len = 0.2,
              bar_thick = 0.030,
              tick_dy = 0.9,    # move numbers closer to bar
              tick_dx_left = -0.18, tick_dx_right = 0.22,
              tick_cex = 1.5,   # ← bigger tick numbers
              title_cex = 1.50,
              title_dy = 1.2)

draw_colorbar(attr(gloriTrack, "pal"),
              expression("GLORI " * m^6 * "A/A (%)"),
              x = 0.5, y = 0.3, orientation = "horizontal",
              bar_len = 0.2, bar_thick = 0.030,
              tick_dy = 0.9,    # move numbers closer to bar
              tick_dx_left = -0.18, tick_dx_right = 0.22,
              tick_cex = 1.5,   # ← bigger tick numbers
              title_cex = 1.50, 
              title_dy = 1.2)

dev.off()



#### add boxplots

try({
  font_add(family = "Calibri",
           regular    = "C:/Windows/Fonts/calibri.ttf",
           bold       = "C:/Windows/Fonts/calibrib.ttf",
           italic     = "C:/Windows/Fonts/calibrii.ttf",
           bolditalic = "C:/Windows/Fonts/calibriz.ttf")
  showtext_auto()
}, silent = TRUE)

dorado <- fread("dorado_genome_csv.txt")
glori  <- fread("glori_genome_csv.txt")

for (dt in list(dorado, glori)) {
  if (max(dt$stoich, na.rm = TRUE) <= 1.5) dt[, stoich := stoich * 100]}

dorado[, method := "Dorado"]
glori[, method := "GLORI"]

# Combine datasets
dt_all <- rbindlist(list(dorado, glori), fill = TRUE)
dt_all <- dt_all[!is.na(stoich) & stoich >= 0 & stoich <= 100]

dt_all[, method := factor(method, levels = unique(method))]

# Color palettes
dorado_pal <- c("#FFF2E0", "#FFDAB3", "#E56A00", "#7A3300")
glori_pal  <- c("#5CB6D3", "#2C89A7", "#0E5A71")

# ---- Plot ----
p <- ggplot() +
  # Dorado points + boxplot (orange)
  geom_boxplot(
    data = dt_all[grepl("Dorado", method)],
    aes(x = method, y = stoich),
    width = 0.4, color = "black", fill = "white",
    size = 0.6, outlier.shape = NA  ) +
  geom_jitter(
    data = dt_all[grepl("Dorado", method)],
    aes(x = method, y = stoich, color = stoich),
    width = 0.1, alpha = 0.9, size = 3  ) +
  scale_color_gradientn(colors = dorado_pal, limits = c(0, 100), name = NULL, guide = "none") +
  
  # GLORI points + boxplot (blue)
  new_scale_color() +  # from 'ggnewscale' package if available; if not, just re-add color scale manually
  geom_boxplot(
    data = dt_all[grepl("GLORI", method)],
    aes(x = method, y = stoich),
    width = 0.4, color = "black", fill = "white",
    size = 0.6, outlier.shape = NA  ) +
  geom_jitter(
    data = dt_all[grepl("GLORI", method)],
    aes(x = method, y = stoich, color = stoich),
    width = 0.1, alpha = 0.9, size = 3  ) +
  scale_color_gradientn(colors = glori_pal, limits = c(0, 100), name = NULL, guide = "none") +
  
  # Axes & labels
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    title = NULL,
    x = NULL,
    y = expression("m"^6*"A/A site ratio (%)")  ) +
  theme_minimal(base_family = "Calibri", base_size = 20) +
  theme(
    plot.title  = element_text(size = 30, hjust = 0.8, color = "black"),
    axis.text   = element_text(size = 25, color = "black"),
    axis.title  = element_text(size = 25, color = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.ticks.y = element_line(color = "black", size = 0.8),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(0.3, "lines")  )

# ---- Save ----
ggsave("SPEN_genomic_boxplot.pdf", p, width = 5, height = 5, useDingbats = FALSE)

