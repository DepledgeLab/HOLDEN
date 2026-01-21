suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicFeatures)
  library(rtracklayer)
  library(GenomicRanges)
  library(data.table)
  library(scales)
  library(grid)
  library(showtext)
})

try({
  font_add(family = "Calibri",
           regular = "C:/Windows/Fonts/calibri.ttf",
           bold    = "C:/Windows/Fonts/calibrib.ttf",
           italic  = "C:/Windows/Fonts/calibrii.ttf",
           bolditalic = "C:/Windows/Fonts/calibriz.ttf")
  showtext_opts(dpi = 300)
  showtext_auto(enable = TRUE)
}, silent = TRUE)

track_font <- "Calibri"

if (!requireNamespace("Cairo", quietly = TRUE)) install.packages("Cairo")

try({
  font_add(family = "Calibri",
           regular = "C:/Windows/Fonts/calibri.ttf",
           bold    = "C:/Windows/Fonts/calibrib.ttf",
           italic  = "C:/Windows/Fonts/calibrii.ttf",
           bolditalic = "C:/Windows/Fonts/calibriz.ttf")
  showtext_auto()
}, silent = TRUE)


gff_file   <- "SPEN.gff3"        
sites_file <- "dorado_tx_csv.txt"  

box_bp     <- 40L       
padding_bp <- 400L    
track_font <- "Calibri" 
label_cex <- track_label_cex    

title_fontsize   <- 36  
track_label_cex  <- 3.0  
track_text_cex   <- 1.6  
axis_cex         <- 1.4  
colorbar_cex     <- 12   
m6A_color <- "#E56A00"

strip_ver <- function(x) sub("\\.\\d+$", "", x)

exons_for_tx <- function(txdb, tx_id) {
  ex_by_tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  if (!tx_id %in% names(ex_by_tx)) {
    id_nover <- strip_ver(tx_id)
    hit <- names(ex_by_tx)[strip_ver(names(ex_by_tx)) == id_nover]
    if (length(hit)) {
      tx_id <- hit[1]
    } else {
      return(NULL)
    }  }
  gr <- ex_by_tx[[tx_id]]
  if (as.character(unique(strand(gr))) == "-") {
    gr <- gr[order(end(gr), decreasing = TRUE)]
  } else {
    gr <- gr[order(start(gr), decreasing = FALSE)]
  }
  gr}

# Map transcript-relative position (1-based) to genomic coordinate (1-based)
map_txpos_to_genome <- function(tx_exons, tx_pos) {
  if (is.null(tx_exons) || length(tx_exons) == 0) return(NULL)
  tx_strand <- as.character(unique(strand(tx_exons)))
  exon_widths <- width(tx_exons)
  cum_end   <- cumsum(exon_widths)
  cum_start <- c(1, head(cum_end, -1) + 1)
  
  hit <- which(tx_pos >= cum_start & tx_pos <= cum_end)
  if (!length(hit)) return(NULL)
  i <- hit[1]
  pos_within_exon <- tx_pos - cum_start[i] + 1
  
  if (tx_strand == "+") {
    g_start <- start(tx_exons)[i] + pos_within_exon - 1L
  } else {
    g_start <- end(tx_exons)[i]   - pos_within_exon + 1L
  }
  list(genomic_start = as.integer(g_start),
       strand        = tx_strand,
       chr           = as.character(seqnames(tx_exons)[i]))}

# Keep boxes compact; we draw start..start+box_bp-1
clamp_box_end <- function(g_start, box_bp, tx_strand) {
  c(g_start, g_start + box_bp - 1L)}

draw_colorbar <- function(pal_fun, title = expression("m"^6*"A/A (%)"),
                          x = 0.80, y = 0.08, width = 0.16, height = 0.03) {
  vp <- viewport(x = x, y = y, width = width, height = height, just = c("left","bottom"))
  pushViewport(vp)
  n  <- 140
  xs <- seq(0, 1, length.out = n)
  for (i in 1:(n-1)) {
    grid.rect(x = (xs[i]+xs[i+1])/2, y = 0.5,
              width = xs[i+1]-xs[i], height = 1,
              gp = gpar(fill = pal_fun(100 * xs[i]), col = NA))  }
  grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.7))
  grid.text(title, x = 0.8, y = 1.35, gp = gpar(fontsize = 10, fontfamily = track_font))
  grid.text("0",   x = 0,   y = -0.25, gp = gpar(fontsize = 9,  fontfamily = track_font))
  grid.text("100", x = 1,   y = -0.25, gp = gpar(fontsize = 9,  fontfamily = track_font))
  popViewport()}


txdb <- makeTxDbFromGFF(gff_file, format = "gff3")

sites <- fread(sites_file, sep = "\t", header = TRUE)
setnames(sites, tolower(names(sites))) 

req <- c("transcript_id", "tx_pos", "strand", "stoich")
stopifnot(all(req %in% names(sites)))
has_name <- "name" %in% names(sites)

tx_order <- unique(sites$transcript_id)

first_non_na <- function(x) x[match(TRUE, !is.na(x))]
if (has_name) {
  tx_label_map <- sites[,
                        .(label = { v <- name; if (all(is.na(v))) unique(transcript_id)[1] else first_non_na(v) }),
                        by = transcript_id  ]
} else {
  tx_label_map <- sites[, .(label = unique(transcript_id)[1]), by = transcript_id]}

site_gr_list <- list()
view_chr <- NULL; view_starts <- c(); view_ends <- c()

for (tx in tx_order) {
  tx_ex <- exons_for_tx(txdb, tx)
  if (is.null(tx_ex)) {
    warning(sprintf("Transcript %s not found in GFF3; skipping.", tx))
    next  }
  
  dtx <- sites[transcript_id == tx]
  if (!nrow(dtx)) next
  
  mapped <- lapply(seq_len(nrow(dtx)), function(i) {
    m <- map_txpos_to_genome(tx_ex, as.integer(dtx$tx_pos[i]))
    if (is.null(m)) return(NULL)
    be <- clamp_box_end(m$genomic_start, box_bp, m$strand)
    data.table(chr   = m$chr,
               start = be[1],
               end   = be[2],
               strand= m$strand,
               stoich= as.numeric(dtx$stoich[i]),
               tx_id = tx)  })
  mapped <- rbindlist(Filter(Negate(is.null), mapped))
  if (!nrow(mapped)) next
  
  site_gr_list[[tx]] <- GRanges(seqnames = mapped$chr,
                                ranges   = IRanges(mapped$start, mapped$end),
                                strand   = mapped$strand,
                                stoich   = mapped$stoich,
                                tx_id    = mapped$tx_id)

  view_chr    <- as.character(seqnames(tx_ex)[1])
  view_starts <- c(view_starts, min(start(tx_ex)))
  view_ends   <- c(view_ends,   max(end(tx_ex)))}

if (!length(site_gr_list)) stop

from <- min(view_starts) - 100
to   <- max(view_ends)   + 100
chr  <- view_chr

overlay_tracks <- list()

for (tx in names(site_gr_list)) {
  tx_ex <- exons_for_tx(txdb, tx)
  tx_label <- tryCatch(tx_label_map[transcript_id == tx, label][1], error = function(e) tx)
  if (is.na(tx_label) || tx_label == "") tx_label <- tx

  geneTrack <- GeneRegionTrack(
    tx_ex,
    chromosome = chr,
    name       = tx_label,
    showId     = FALSE,
    collapseTranscripts = TRUE,     # ensures single collapsed transcript view
    transcriptAnnotation = "id",
    cex        = 1.0,
    fontfamily = track_font,
    fill       = "#EAEAEA",        
    col        = "#444444",        
    lwd        = 0.8,
    col.line   = "#444444",       
    col.arrow  = "#444444",       
    lty        = 1,                
    arrowHeadWidth = 6,           
    thinBoxFeature = TRUE   )

  sg <- site_gr_list[[tx]]
  
  siteTrack <- AnnotationTrack(
    range = sg,
    chromosome = chr,
    name = "",
    shape = "box",
    stacking = "dense",
    fill = m6A_color,    
    col  = m6A_color,  
    cex  = 0.6,
    fontfamily = track_font  )

  otrack <- OverlayTrack(trackList = list(geneTrack, siteTrack), name = tx_label)
  displayPars(otrack) <- list(
    rotation.title = 0,
    cex.title      = label_cex,
    fontfamily     = track_font  )
  
  overlay_tracks[[length(overlay_tracks) + 1L]] <- otrack}

Cairo::CairoPDF("SPEN_isoforms_m6A_boxes.pdf", width = 20, height = 6)

grid::grid.newpage()

title_offset_top   <- 1.0  
title_to_tracks_gap <- 2.2  

grid::grid.text(
  expression(italic("SPEN") ~"isoforms with " * m^6 * "A sites" ~(Dorado)),
  x = 0.58,
  y = grid::unit(1, "npc") - grid::unit(0.6, "lines"),
  just = "center",
  gp = grid::gpar(fontsize = title_fontsize,
                  fontfamily = track_font,
                  fontface = "bold"))

grid::pushViewport(grid::viewport(
  x = grid::unit(0.06, "npc"),     
  width = grid::unit(0.94, "npc"),
  y = grid::unit(0.00, "npc"),
  height = grid::unit(0.88, "npc"),
  just = c("left", "bottom")))

plotTracks(
  overlay_tracks,
  chromosome = chr, from = from, to = to,
  add = TRUE,
  just.title   = "right",
  title.width  = 4.5,
  title.offset = 4,
  background.title = "#FFFFFF",
  col.title    = "#000000",
  cex.title    = track_label_cex,  # enlarged isoform names
  fontfamily   = track_font,
  rotation.title = 0,
  min.distance = 0.1,
  cex.axis     = axis_cex,         # NEW: axis label size
  cex          = track_text_cex,   # NEW: text inside tracks
  sizes        = rep(2, length(overlay_tracks)))

grid::popViewport()

dev.off()



# boxplot: stoichiometry for SPEN-201 vs SPEN-202 ---


library(data.table)
library(ggplot2)
library(showtext)
library(scales)

try({
  font_add(family = "Calibri",
           regular    = "C:/Windows/Fonts/calibri.ttf",
           bold       = "C:/Windows/Fonts/calibrib.ttf",
           italic     = "C:/Windows/Fonts/calibrii.ttf",
           bolditalic = "C:/Windows/Fonts/calibriz.ttf")
  showtext_auto()
}, silent = TRUE)

dt <- fread("dorado_tx_csv.txt")

dt_sub <- dt[name %in% c("SPEN-201", "SPEN-202")]


# ---- Color scale for stoichiometry (same as Gviz gradient) ----
stoich_pal <- c("#FFF2E0", "#FFDAB3", "#E56A00", "#7A3300")
col_fun <- scales::col_numeric(stoich_pal, domain = c(0, 100))

# ---- Plot ----
p <- ggplot(dt_sub, aes(x = name, y = stoich)) +
  geom_boxplot(
    width = 0.4,
    color = "black",
    fill  = "white",      # white background for boxes
    size  = 0.5,
    outlier.shape = NA  ) +
  geom_jitter(
    aes(color = stoich),  # color dots by stoichiometry
    width = 0.1,
    alpha = 0.9,
    size = 3
  ) +
  scale_color_gradientn(
    colors = stoich_pal,
    limits = c(0, 100)
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = expression("m"^6*"A/A site ratio (%)")
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.05))) +
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

ggsave("SPEN_isoforms_boxplot.pdf", p, width = 5, height = 5, useDingbats = FALSE)


