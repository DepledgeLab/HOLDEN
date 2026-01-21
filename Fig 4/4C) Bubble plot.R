## DAVID analysis was conducted using 20 genes minimum 

library(data.table)
library(ggplot2)
library(data.table)

overlap <- fread("DAVID_GO_NHDFs_overlap_HD106.txt")
NHDF    <- fread("DAVID_GO_NHDFs_unique_sites.txt")
HD106   <- fread("DAVID_GO_HD106_unique_sites.txt")

FDR_MAX <- 0.01
TOP_N   <- 20
dataset_order <- c("overlap", "NHDF", "HD106")

coerce_numeric <- function(dt) {
  num_cols <- c("Count","PValue","Fold_Enrichment","Bonferroni","Benjamini","FDR",
                "%","List Total","Pop Hits","Pop Total")
  for (cn in intersect(names(dt), num_cols)) {
    suppressWarnings(dt[[cn]] <- as.numeric(dt[[cn]]))
  }  dt}

overlap <- coerce_numeric(copy(overlap))
NHDF    <- coerce_numeric(copy(NHDF))
HD106   <- coerce_numeric(copy(HD106))

ov <- copy(overlap)
ov[, Term := sub("^.*~", "", Term)]
ov_filt <- ov[is.finite(FDR) & FDR <= FDR_MAX]

setorder(ov_filt, -Count, -`Fold_Enrichment`)
ov_top <- ov_filt[1:min(TOP_N, .N)]

term_levels <- ov_top[order(-`Count`), unique(Term)]

prep_ds <- function(dt, label) {
  dt <- copy(dt)
  dt[, Term := sub("^.*~", "", Term)]
  dt <- dt[Term %in% term_levels]
  dt[, Dataset := label]
  dt[, log10FDR := -log10(pmax(FDR, 1e-300))]
  dt}

plot_df <- rbindlist(list(
  prep_ds(ov_top, "Shared"), 
  prep_ds(NHDF,   "NHDF"),
  prep_ds(HD106,  "HD106")
), use.names = TRUE, fill = TRUE)

dataset_order <- c("Shared", "NHDF", "HD106")

plot_df[, Term := factor(Term, levels = rev(term_levels))]
plot_df[, Dataset := factor(Dataset, levels = dataset_order,
                            labels = c("Shared", "NHDF", "HD10.6"))]

p <- ggplot(plot_df, aes(x = Dataset, y = Term)) +
  geom_point(aes(size = Count, color = log10FDR), alpha = 0.9) +
  scale_size_area(max_size = 12, name = "Count") +
  scale_color_gradient(
    name = expression(-log[10](FDR)),
    low  = "#3B4CC0",  # blue
    high = "#D63C4A"   # red
  ) +
  labs(
    x = NULL, y = NULL,
    title = "GO term analysis (Biological process)",
    subtitle = sprintf("FDR ≤ %.2f, top %d by Count", 
                       FDR_MAX, length(term_levels))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y   = element_text(size = 34),
    axis.text.x   = element_text(size = 34, face = "bold"),
    legend.title  = element_text(size = 30, face = "bold"),
    legend.text   = element_text(size = 30),
    plot.title    = element_text(size = 32, face = "bold"),
    plot.subtitle = element_text(size = 30),
    plot.margin   = margin(14, 10, 14, 10)  )

ggsave(
  filename = "GO_bubble_overlap_top20.pdf",
  plot     = p,
  width    = 23, height = 13, units = "in",
  device   = grDevices::cairo_pdf)
