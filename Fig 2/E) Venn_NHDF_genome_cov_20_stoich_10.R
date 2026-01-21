#install.packages("data.table")
#install.packages("readr")
#install.packages("VennDiagram")
#install.packages("readr")
#install.packages("extrafont")

library(data.table)
library(VennDiagram)
library(grid)
library(extrafont)

#extrafont::loadfonts(device = "pdf")  # once per session
#loadfonts(device = "pdf")  # Use "pdf" for saving PDF outputs
#fonts()  # This will display a list of all available fonts

# Load all four tables
table1 <- fread("a.NHDF_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.Filt.motif.20reads_0stoich.bed")
table2 <- fread("a.NHDF_DMSO_48h_1.sup-allMods.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.inosine.noFilt.motif.20reads_0stoich.bed")
table3 <- fread("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.Filt.motif.20reads_0stoich.bed")
table4 <- fread("NHDF_DMSO_48h_1.sup-m6A_DRACH.trimAdapters.dorado.0.9.0.hg38.sorted.m6A.noFilt.motif.20read_0stoich.bed")

table1 <- table1[V11 >= 10]
table2 <- table2[V11 >= 10]
table3 <- table3[V11 >= 10]
table4 <- table4[V11 >= 10]

setnames(table1, old = c("V1", "V2"), new = c("chrom", "start_position"))
setnames(table2, old = c("V1", "V2"), new = c("chrom", "start_position"))
setnames(table3, old = c("V1", "V2"), new = c("chrom", "start_position"))
setnames(table4, old = c("V1", "V2"), new = c("chrom", "start_position"))

table1$pair <- paste(table1$chrom, table1$start_position, sep=":")
table2$pair <- paste(table2$chrom, table2$start_position, sep=":")
table3$pair <- paste(table3$chrom, table3$start_position, sep=":")
table4$pair <- paste(table4$chrom, table4$start_position, sep=":")

venn_list <- list(
  table1 = unique(table1$pair),
  table2 = unique(table2$pair),
  table3 = unique(table3$pair),
  table4 = unique(table4$pair))

pdf("Venn_NHDF_genome.pdf", width = 12, height = 12, family = "Calibri", useDingbats = FALSE)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,  # Do not save to file at this point
  output = TRUE,    # Output to current device
  lwd = 1,
  fill = c("#FFDAB9", "gray", "#FFB080", "dimgray"),
  alpha = 0.5,
  cat.cex = 2.3,
  cex = 2.8,
  margin   = 0.06,
  
  category.names = c(
    expression(atop("all-context model", "filtered")),     # table1
    expression(atop("all-context model", "unfiltered")),   # table2
    expression(atop("DRACH model", "filtered")),           # table3
    expression(atop("DRACH model", "unfiltered"))          # table4
  ))

for (j in seq_along(venn.plot)) {
  if (inherits(venn.plot[[j]], "text")) {
    lab <- venn.plot[[j]]$label
    if (grepl("^[0-9]+$", lab)) {
      venn.plot[[j]]$label <- format(as.numeric(lab), big.mark = ",", scientific = FALSE)
    }  }}

grid.draw(venn.plot)
dev.off()


