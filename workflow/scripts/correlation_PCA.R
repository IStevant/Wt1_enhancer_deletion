source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("grid")
  library("viridis")
  library("ggplot2")
  library("ComplexHeatmap")
  library("dplyr")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

norm_data <- read.csv(file = snakemake@input[["norm_counts"]], header=TRUE, row.names = 1, check.names = FALSE)

conditions <- gsub("_Rep.*", "\\1", colnames(norm_data))

# Order the samples in the graphs according to their order in the expression matrix
cond_sorted <- unique(conditions)

# Pick colors
set.seed(1234)
# colours <- randomcoloR::distinctColorPalette(length(cond_sorted))
colours <- c(
  "#2f8e26",
  "#83da7a",
  "#d5f3d2",
  "#9542ca",
  "#cba2e5",
  "#ebdcf6"
)
names(colours) <- cond_sorted

corr_method <- snakemake@params[["corr_method"]]
cluster_method <- "complete"

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the correlation matrix between the samples
#' @param matrix Expression matrix.
#' @param method Method for the correlation (Example: "Spearman" or "Pearson"). Default is "Spearman".
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @return Pheatmap object.
correlation <- function(matrix, method = "Spearman", colours, cluster_method="ward.D2") {
  # matrix <- matrix[, order(names(matrix))]
  cor_data <- cor(matrix, method = method)
  cor_data[cor_data == 1.000] <- NA
  conditions <- factor(conditions, levels = cond_sorted)

  anno <- data.frame(
    Samples = conditions
  )

  colors <- list(Samples=colours)

  heatmap <- Heatmap(
    cor_data,
    name = "Correlation",
    # cluster_rows = F, 
    # cluster_columns = F, 
    clustering_method_columns = cluster_method,
    clustering_method_rows = cluster_method,
    # column_split = 3,
    # row_split = 3,
    column_title = NULL,
    row_title = NULL,
    col = viridis::viridis(n = 100,option = 'C'),
    left_annotation = rowAnnotation(df = anno, col=colors),
    top_annotation = columnAnnotation(df = anno, col=colors, show_legend=c(F,F,F)),
    width = unit(0.65, "snpc"),
    height = unit(0.65, "snpc"),
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14)
  )

  heatmap <- draw(heatmap,
      column_title=paste0("Pairwise sample correlation (",method,")"),
      column_title_gp=grid::gpar(fontsize=14, fontface="bold")
    )

  plot <- grid.grabExpr(
    draw(heatmap)
  )

  return(plot)
}

#' Proceed to the PCA using prcomp.
#' @param matrix Expression matrix.
#' @return prcomp object.
run_pca <- function(matrix) {
  print("Calculating the PCA...")
  t.matrix <- t(matrix)
  t.matrix.no0 <- t.matrix[, colSums(t.matrix) != 0]
  pca <- prcomp(
    t.matrix.no0,
    center = TRUE,
    scale. = TRUE
  )
  return(pca)
}

#' Plot the PCA.
#' @param matrix Expression matrix.
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @param PCs Vector PCs to plot. Default is "PC1" vs "PC2". If more than two PCs provided, all the possible combinations are drawn.
#' @return Pheatmap object.
plot.pca <- function(matrix, colours, PCs = c("PC1", "PC2")) {
  pca <- run_pca(matrix)
  print("Generating the plots...")
  percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  cond <- factor(conditions)
  col <- factor(conditions)
  levels(col) <- colours
  col <- as.vector(col)
  replicates <- stringr::str_to_sentence(sub(".*_(Rep[0-9]+)", "\\1", colnames(matrix)))
  sex <- sapply(strsplit(colnames(matrix), "_"), `[`, 1)
  scores <- as.data.frame(pca$x)
  PCs.combinations <- combn(PCs, 2)
  plots <- apply(
    PCs.combinations,
    2,
    function(combination) {
      data <- as.data.frame(scores[, c(combination[1], combination[2])])
      data$cond <- conditions
      data$rep <- replicates
      colnames(data) <- c("PC_x", "PC_y", "cond", "rep")
      data$cond = factor(data$cond, levels = cond_sorted)

      plot <- ggplot(data, aes(x = PC_x, y = PC_y, fill = cond, color = cond)) +
        geom_point(aes(shape = sex), size = 6) +
        scale_fill_manual(values = colours) +
        scale_color_manual(values = colours) +
        # ggtitle("PCA on TPM")+
        theme_bw() +
        xlab(paste(combination[1], " ", "(", round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit = 2), "%)", sep = "")) +
        ylab(paste(combination[2], " ", "(", round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit = 2), "%)", sep = "")) +
        guides(
          fill = guide_legend(override.aes = list(shape = 21)),
          shape = guide_legend(override.aes = list(fill = "#444444"))
        ) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          aspect.ratio = 1,
          legend.title = element_blank(),
          legend.text = element_text(size = 14)
        )
      return(plot)
    }
  )
  print("Done.")
  return(plots)
}

#################################################################################################################################

###########################################
#                                         #
#        Plot correlation and PCA         #
#                                         #
###########################################


corr_plot <- correlation(norm_data, corr_method, colours, cluster_method)
pca_plot <- plot.pca(norm_data, colours, c("PC1", "PC2"))

# Combine the two plots
figure <- plot_grid(
  plotlist = list(corr_plot, pca_plot[[1]]),
  labels = "AUTO",
  ncol = 2, 
  align = "h"
)

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 55,
  base_height = 20,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 55,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

# Remove the unwanted file
file.remove("Rplots.pdf")