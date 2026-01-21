library(data.table)
library(VennDiagram)
library(grid)
library(sysfonts)
library(showtext)
library(ggplot2)

font_add("Calibri",
         regular = "C:/Windows/Fonts/calibri.ttf",
         bold    = "C:/Windows/Fonts/calibrib.ttf")
showtext_auto()

table1 <- fread("hg38_NHDF_DRACH_highConfidence_with_gene.txt")
table2 <- fread("hg38_HD106_DRACH_highConfidence_with_gene.txt")
table2[, chrom := sub("^chr", "", chrom)]

table1 <- table1[!is.na(gene) & gene != ""]
table2 <- table2[!is.na(gene) & gene != ""]

common_genes <- intersect(unique(table1$gene), unique(table2$gene))
table1 <- table1[gene %in% common_genes]
table2 <- table2[gene %in% common_genes]

table1[, pair := paste(chr, start, sep=":")]
table2[, pair := paste(chrom, start_position, sep=":")]

venn_list <- list(
  NHDF  = unique(table1$pair),
  HD10.6 = unique(table2$pair))


venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,    
  lwd = 1.5,
  fill = c("#FFDAB9", "#A0D3E8"),
  alpha = 0.55,
  cex = 0,
  cat.cex = 3.5,           
  fontfamily = "Calibri",
  cat.fontfamily = "Calibri",
  cat.pos = c(-45, 45),    
  cat.dist = c(0.05, 0.05),
  margin = 0.1)

for (j in seq_along(venn.plot)) {
  if (inherits(venn.plot[[j]], "text")) {
    lab <- venn.plot[[j]]$label
    if (grepl("^[0-9]+$", lab)) {
      venn.plot[[j]]$label <- format(as.numeric(lab), big.mark = ",", scientific = FALSE)
    }}}

pdf("Venn_dorado_cellTypes.pdf", width = 10, height = 10, useDingbats = FALSE)
grid.newpage()
grid.draw(venn.plot)

grid.text(
  expression("High-confidence " * m^6 * "A sites"),
  x = 0.5, y = 0.85,
  gp = gpar(fontfamily = "Calibri", fontsize = 45))

dev.off()

