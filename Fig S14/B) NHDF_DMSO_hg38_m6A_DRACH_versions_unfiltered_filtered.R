library(UpSetR)

#Unfiltered
m6A_DRACH_Unfiltered <- c(
  Unfiltered_0.5.0 = 1204,
  Unfiltered_0.6.0 = 6915,
  Unfiltered_0.9.0 = 1,
  "Unfiltered_0.5.0&Unfiltered_0.6.0" = 17777,
  "Unfiltered_0.8.0&Unfiltered_0.9.0" = 7453,
  "Unfiltered_0.5.0&Unfiltered_0.8.0&Unfiltered_0.9.0" = 365,
  "Unfiltered_0.6.0&Unfiltered_0.8.0&Unfiltered_0.9.0" = 1154,
  "Unfiltered_0.5.0&Unfiltered_0.6.0&Unfiltered_0.8.0&Unfiltered_0.9.0" = 122679
)


# Plot
# empty.intersections = "on", 
DRACH_un <-upset(fromExpression(m6A_DRACH_Unfiltered), order.by = "freq", matrix.color="red", 
                 decreasing = T, mb.ratio = c(0.65, 0.35), text.scale = 2, point.size = 3, 
                 line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "red", 
                 main.bar.color = "red"
)
DRACH_un

pdf(file = "B) NHDF_DMSO_hg38_m6A_DRACH_versions_unfiltered_upsetplot.pdf", width = 8, height = 5)
DRACH_un
dev.off()


#Filtered
m6A_DRACH_filtered <- c(
  Filtered_0.5.0 = 1084,
  Filtered_0.6.0 = 653,
  Filtered_0.9.0 = 1,
  "Filtered_0.5.0&Filtered_0.6.0" = 6687,
  "Filtered_0.8.0&Filtered_0.9.0" = 7711,
  "Filtered_0.5.0&Filtered_0.8.0&Filtered_0.9.0" = 545,
  "Filtered_0.6.0&Filtered_0.8.0&Filtered_0.9.0" = 432,
  "Filtered_0.5.0&Filtered_0.6.0&Filtered_0.8.0&Filtered_0.9.0" = 109685
)


# Plot
# empty.intersections = "on", 
DRACH_filt <-upset(fromExpression(m6A_DRACH_filtered), order.by = "freq", matrix.color="red", 
                   decreasing = T, mb.ratio = c(0.65, 0.35), text.scale = 2, point.size = 3, 
                   line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "red", 
                   main.bar.color = "red"
)
DRACH_filt



pdf(file = "B) NHDF_DMSO_hg38_m6A_DRACH_versions_filtered_upsetplot.pdf", width = 8, height = 5)
DRACH_filt
dev.off()
