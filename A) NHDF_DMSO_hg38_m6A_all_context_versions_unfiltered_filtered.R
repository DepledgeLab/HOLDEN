library(UpSetR)

#Unfiltered
m6A_all_context_Unfiltered <- c(
  Unfiltered_0.7.0 = 397863,
  Unfiltered_0.8.0 = 3,
  Unfiltered_0.9.0 = 12,
  "Unfiltered_0.7.0&Unfiltered_0.8.0" = 1,
  "Unfiltered_0.7.0&Unfiltered_0.9.0" = 5,
  "Unfiltered_0.8.0&Unfiltered_0.9.0" = 70725,
  "Unfiltered_0.7.0&Unfiltered_0.8.0&Unfiltered_0.9.0" = 373903
)


# Plot
# empty.intersections = "on", 
all_context_un <- upset(fromExpression(m6A_all_context_Unfiltered), order.by = "freq", matrix.color="red", 
                        decreasing = T, mb.ratio = c(0.7, 0.3), text.scale = 2, point.size = 3, 
                        line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "red", 
                        main.bar.color = "red"
)
all_context_un

pdf(file = "A) NHDF_DMSO_hg38_m6A_all_context_0.7.0_0.8.0_0.9.0_unfiltered_upsetplot.pdf", width = 8, height = 5)
all_context_un
dev.off()



#Filtered
m6A_all_context_Filtered <- c(
  Filtered_0.7.0 = 16752,
  Filtered_0.8.0 = 0,
  Filtered_0.9.0 = 3,
  "Filtered_0.7.0&Filtered_0.8.0" = 0,
  "Filtered_0.7.0&Filtered_0.9.0" = 0,
  "Filtered_0.8.0&Filtered_0.9.0" = 7407,
  "Filtered_0.7.0&Filtered_0.8.0&Filtered_0.9.0" = 127291
)

# Plot
# empty.intersections = "on", 
all_context_filt <-upset(fromExpression(m6A_all_context_Filtered), order.by = "freq", matrix.color="red", 
                         decreasing = T, mb.ratio = c(0.7, 0.3), text.scale = 2, point.size = 3, 
                         line.size = 1.5, show.numbers = FALSE, sets.x.label = "", sets.bar.color = "red", 
                         main.bar.color = "red"
)
all_context_filt


pdf(file = "A) NHDF_DMSO_hg38_m6A_all_context_0.7.0_0.8.0_0.9.0_filtered_upsetplot.pdf", width = 8, height = 5)
all_context_filt
dev.off()
