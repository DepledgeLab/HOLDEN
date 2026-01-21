library(data.table)
library(ggplot2)
library(showtext)

font_add(family = "Calibri",
         regular = "C:/Windows/Fonts/calibri.ttf",
         bold    = "C:/Windows/Fonts/calibrib.ttf",
         italic  = "C:/Windows/Fonts/calibrii.ttf",
         bolditalic = "C:/Windows/Fonts/calibriz.ttf")
showtext_auto()

dt <- fread("DAVID_KEGG_NHDFs_overlap_HD106.txt", sep = "auto", header = TRUE, quote = "")

plot_df <- dt[, .(
  Term = sub("^.*~", "", as.character(Term)),
  Fold_Enrichment = as.numeric(Fold_Enrichment),
  FDR = as.numeric(FDR),
  Count = as.numeric(Count))]

plot_df <- plot_df[is.finite(Fold_Enrichment) & is.finite(FDR) & is.finite(Count) & !is.na(Term)]
plot_df <- plot_df[FDR <= 0.01]

plot_df <- plot_df[order(-Count, -Fold_Enrichment)][1:min(20, .N)]

plot_df <- plot_df[order(-Count, -Count)]
plot_df[, Term := factor(Term, levels = rev(Term))]

plot_df[, NegLog10FDR := -log10(pmax(FDR, 1e-300))]

p <- ggplot(plot_df, aes(Fold_Enrichment, Term)) +
  geom_segment(aes(x = 0, xend = Fold_Enrichment, y = Term, yend = Term),
               linewidth = 0.6, color = "grey60") +
  geom_point(aes(size = Count, color = NegLog10FDR)) +
  scale_size(range = c(2.5, 8), name = "Count") +
  scale_color_gradient(name = expression(-log[10](FDR)),
                       low = "#3B4CC0", high = "#D63C4A") +
  labs(x = "Fold enrichment", y = NULL,
       title = "KEGG pathway analysis",
       subtitle = expression("38,217 shared" ~ m^6 * "A sites (7,019 genes)")
  ) +
  
  theme_minimal(base_size = 12, base_family = "Calibri") +
  theme(
    axis.text.y     = element_text(size = 28, family = "Calibri"),
    axis.text.x     = element_text(size = 28, family = "Calibri"),
    axis.title.x    = element_text(size = 28, family = "Calibri"),
    legend.title    = element_text(size = 24, family = "Calibri"),
    legend.text     = element_text(size = 24, family = "Calibri"),
    plot.title    = element_text(size = 30, face = "bold", family = "Calibri"),
    plot.subtitle = element_text(size = 24, family = "Calibri"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
ggsave("KEGG_analysis_overlapping_genes_top20.pdf", p, width = 14, height = 11, useDingbats = FALSE)

## DAVID KEGG analysis was conducted using 20 genes minimum

### correlation plots for single genes (example WNTB5)

dt <- fread("NHDF_HD106_overlap.txt", sep = "auto", header = TRUE, quote = "")

WNT5B_dt <- dt[grepl("WNT5B", `gene.NHDF`, fixed = TRUE)]

p <- ggplot(WNT5B_dt, aes(x = `percent_modified_DMSO-NHDF`, 
                          y = `percent_modified_DMSO-HD106`)) +
  geom_point(size = 5, alpha = 1, color = "darkblue") +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "solid", size = 2, color = "black") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    title = expression(bolditalic("WNT5B")~ bold("shared " * m^6 * "A sites")),
    x = expression(m^6 * "A/A site ratio NHDF (%)"),
    y = expression(m^6 * "A/A site ratio HD10.6 (%)")
  ) +
  theme_minimal(base_size = 9, base_family = "Calibri") +
  theme(
    plot.title  = element_text(size = 39, face = "bold", hjust = 0.5, family = "Calibri", margin = margin(b = 20) ),
    axis.title  = element_text(size = 38, family = "Calibri"),
    axis.text   = element_text(size = 38, family = "Calibri"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(t = 20, r = 30, b = 60, l = 50)
  )

ggsave("WNT5B_scatter.pdf", p, width = 7.5, height = 7.5, units = "in", dpi = 300)
