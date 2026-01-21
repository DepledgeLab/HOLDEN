library(UpSetR)

#Unfiltered
pseU_Unfiltered <- c(
  Unfiltered_0.7.0 = 234244,
  Unfiltered_0.8.0 = 6,
  Unfiltered_0.9.0 = 4,
  "Unfiltered_0.7.0&Unfiltered_0.8.0" = 0,
  "Unfiltered_0.7.0&Unfiltered_0.9.0" = 4,
  "Unfiltered_0.8.0&Unfiltered_0.9.0" = 239652,
  "Unfiltered_0.7.0&Unfiltered_0.8.0&Unfiltered_0.9.0" = 228956
)


# Plot
# empty.intersections = "on", 
pseU_un <-upset(fromExpression(pseU_Unfiltered), order.by = "freq", matrix.color="blue", 
                decreasing = T, mb.ratio = c(0.7, 0.3), text.scale = 2, point.size = 3, 
                line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "blue", 
                main.bar.color = "blue"
)
pseU_un


pdf(file = "C) NHDF_DMSO_hg38_pseU_versions_unfiltered_upsetplot.pdf", width = 8, height = 5)
pseU_un
dev.off()


#Filtered
pseU_Filtered <- c(
  Filtered_0.7.0 = 6431,
  Filtered_0.8.0 = 0,
  Filtered_0.9.0 = 0,
  "Filtered_0.7.0&Filtered_0.9.0" = 0,
  "Filtered_0.8.0&Filtered_0.9.0" = 8446,
  "Filtered_0.7.0&Filtered_0.8.0&Filtered_0.9.0" = 3046
)

# Plot
# empty.intersections = "on", 
pseU_filt <-upset(fromExpression(pseU_Filtered), order.by = "freq", matrix.color="blue", 
                  decreasing = T, mb.ratio = c(0.7, 0.3), text.scale = 2, point.size = 3, 
                  line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "blue", 
                  main.bar.color = "blue"
)
pseU_filt



pdf(file = "C) NHDF_DMSO_hg38_pseU_versions_filtered_upsetplot.pdf", width = 8, height = 5)
pseU_filt
dev.off()

