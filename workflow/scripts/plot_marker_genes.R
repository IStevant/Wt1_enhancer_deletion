source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["TPM"]], row.names = 1, check.names = FALSE)
markerGenes <- read.csv(file = snakemake@input[["marker_genes"]])

conditions <- gsub("_Rep.*", "\\1", colnames(TPM))

# Order the samples in the graphs according to their order in the expression matrix
cond_sorted <- unique(conditions)

# Pick colors
set.seed(1234)
conditions_color <- randomcoloR::distinctColorPalette(length(cond_sorted))
names(conditions_color) <- cond_sorted

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the expression plots marker genes for each cell type
#' @param TPM TPM expression matrix.
#' @param markers Dataframe containing the marker genes and their corresponding cell types.
#' @return Ggplot object.
plot_expression_scatter <- function(TPM, markers) {
  cellTypes <- unique(markers$Cell_type)
  plot_list <- lapply(cellTypes, function(cellType) {
    genes <- markers[markers$Cell_type %in% cellType, "Gene"]
    gene_exp(genes, TPM, cellType)
  })
}

#' Draw the expression plots for each marker gene of a given cell type
#' @param genes Vector containing the marker gene names.
#' @param TPM TPM expression matrix.
#' @param title Cell type name.
#' @return Ggplot object.
gene_exp <- function(genes, TPM, title) {
  plotlist <- list()
  for (gene in genes) {
    exp <- as.numeric(TPM[gene, ])
    gene_exp <- data.frame(
      condition = conditions,
      exp = exp
    )
    df <- dplyr::group_by(gene_exp, condition)
    options(dplyr.summarise.inform = FALSE)
    df.summary2 <- dplyr::summarise(
      df,
      sd = sd(exp),
      len = mean(exp)
    )

    df.summary2$condition = factor(df.summary2$condition, levels = cond_sorted)
    plot <- ggplot(df.summary2, aes(x = condition, y = len, fill=condition, color=condition)) +
      geom_bar(stat="identity") +
      geom_jitter(data=gene_exp, aes(x = condition, y = exp), width = 0.25, size = 0.7, color = "#666666") +
      # geom_line(linewidth = 1) +
      # geom_point(size = 2.5) +
      geom_errorbar(
        aes(
          ymin = len - sd,
          ymax = len + sd
        ),
        width = .1,
        color = "black"
      ) +
      scale_fill_manual(
        values = alpha(conditions_color, 0.8)
      ) +
      scale_colour_manual(
        values = conditions_color
      ) +
      scale_y_continuous(labels = scales::comma) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(title = gene, x = "Samples", y = "Expression (TPM)") +
      theme_light() +
      theme(
        plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
        # axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        # axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_blank()
        # aspect.ratio = 0.5
      )
    legend <- get_legend(
      plot + theme(legend.direction = "horizontal", legend.box.margin = margin(0, 0, 0, 0))
    )
    plot <- plot + theme(legend.position = "none")
    plotlist[[gene]] <- plot
  }

  gtitle <- ggdraw() +
    draw_label(
      title,
      fontface = "bold",
      size = 14,
      hjust = 0,
      vjust = 0.5,
      x = 0
    ) +
    theme(
      plot.background = element_rect(fill = "#ECF0F5", color = NA),
      plot.margin = margin(0, 0, 0, 25)
    )
  y_axis <- ggdraw() +
    draw_label(
      "Expression (TPM) vs embryonic stages",
      size = 12,
      hjust = 0.5
    )

  plot_graphs <- plot_grid(plotlist = plotlist, ncol = 4, align = "hv")

  plots <- plot_grid(
    gtitle, plot_graphs,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  return(plots)
}

##########################################
#                                        #
#          Plot gene expression          #
#                                        #
##########################################

# Plot expression of the marker genes
plot_genes <- plot_expression_scatter(TPM, markerGenes)

figure <- plot_grid(
  plotlist = plot_genes,
  labels = "AUTO",
  ncol = 1,
  rel_heights = c(1, 2, 2, 1, 1),
  align = "v"
)


##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["fig_pdf"]],
  figure,
  base_width = 42,
  base_height = 60,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["fig_png"]],
  figure,
  base_width = 42,
  base_height = 60,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

# Remove the unwanted file
file.remove("Rplots.pdf")