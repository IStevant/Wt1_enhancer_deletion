source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
  library("ggpubr")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

TPM <- read.csv(file = snakemake@input[["TPM"]], row.names = 1, check.names = FALSE)
markerGenes <- read.csv(file = snakemake@input[["marker_genes"]])

stat <- snakemake@params[["stat"]]
print("####")
print(stat)
print("####")

conditions <- gsub("_Rep.*", "\\1", colnames(TPM))

# Order the samples in the graphs according to their order in the expression matrix
cond_sorted <- unique(conditions)

# Pick colors
set.seed(1234)
# conditions_color <- randomcoloR::distinctColorPalette(length(cond_sorted))
conditions_color <- c(
  "#2f8e26",
  "#83da7a",
  "#d5f3d2",
  "#9542ca",
  "#cba2e5",
  "#ebdcf6"
)
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
  if (stat == "yes") {

    plotlist <- list()

    for (gene in genes) {
      exp <- as.numeric(TPM[gene, ])
      gene_exp <- data.frame(
        condition = conditions,
        exp = exp
      )

      # Add sex information to gene_exp
      gene_exp$sex <- ifelse(grepl("XX", gene_exp$condition), "XX", ifelse(grepl("XY", gene_exp$condition), "XY", NA))

      # Summarize data
      df <- dplyr::group_by(gene_exp, condition)
      options(dplyr.summarise.inform = FALSE)
      df.summary2 <- dplyr::summarise(
        df,
        sd = sd(exp),
        len = mean(exp)
      )

      # Add sex information to df.summary2
      df.summary2$condition <- factor(df.summary2$condition, levels = cond_sorted)
      df.summary2$sex <- ifelse(grepl("XX", df.summary2$condition), "XX", ifelse(grepl("XY", df.summary2$condition), "XY", NA))

      # Subset gene_exp by sex
      gene_exp_XY <- gene_exp[gene_exp$sex == "XY", ]
      gene_exp_XX <- gene_exp[gene_exp$sex == "XX", ]

      # Perform Wilcoxon tests for XY
      p_XY_Het <- tryCatch({
        t.test(exp ~ condition, data = subset(gene_exp_XY, condition %in% c("XY_WT", "XY_Het")))$p.value
      }, error = function(e) NA)

      p_XY_Hom <- tryCatch({
        t.test(exp ~ condition, data = subset(gene_exp_XY, condition %in% c("XY_WT", "XY_Hom")))$p.value
      }, error = function(e) NA)

      # Perform Wilcoxon tests for XX
      p_XX_Het <- tryCatch({
        t.test(exp ~ condition, data = subset(gene_exp_XX, condition %in% c("XX_WT", "XX_Het")))$p.value
      }, error = function(e) NA)

      p_XX_Hom <- tryCatch({
        t.test(exp ~ condition, data = subset(gene_exp_XX, condition %in% c("XX_WT", "XX_Hom")))$p.value
      }, error = function(e) NA)

      # Create the plot
      plot <- ggplot(df.summary2, aes(x = condition, y = len, fill = condition, color = condition)) +
        geom_bar(stat = "identity") +
        geom_jitter(data = gene_exp, aes(x = condition, y = exp), width = 0.25, size = 0.7, color = "#666666") +
        geom_errorbar(
          aes(
            ymin = len - sd,
            ymax = len + sd
          ),
          width = 0.1,
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
          # axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          aspect.ratio = 0.8
        )

      # Function to add p-value annotations
      add_pvalue <- function(plot, x_start, x_end, y, p_value) {
        if (!is.na(p_value)) {
          plot +
            annotate(
              "segment",
              x = x_start, xend = x_end,
              y = y + 0.1 * max(gene_exp$exp), yend = y + 0.1 * max(gene_exp$exp),
              color = "black", size = 0.5
            ) +
            annotate(
              "text",
              x = (x_start + x_end) / 2,
              y = y + 0.05 * max(gene_exp$exp),
              label = ifelse(p_value < 0.05, paste0("p = ", format.pval(p_value, digits = 2)), "ns"),
              vjust = -1, size = 3.5
            ) +
            annotate(
              "text",
              x = (x_start + x_end) / 2,
              y = y + 0.2 * max(gene_exp$exp),
              label = "",
              vjust = -1, size = 3.5
            )
        } else {
          plot
        }
      }

      # Get x positions for conditions
      x_positions <- as.numeric(factor(df.summary2$condition, levels = levels(df.summary2$condition)))

      # Add p-values for XY comparisons
      if (!is.na(p_XY_Het)) {
        x_start <- x_positions[which(df.summary2$condition == "XY_WT")]
        x_end <- x_positions[which(df.summary2$condition == "XY_Het")]
        plot <- add_pvalue(plot, x_start, x_end, max(gene_exp$exp) + 0.05 * max(gene_exp$exp), p_XY_Het)
      }

      if (!is.na(p_XY_Hom)) {
        x_start <- x_positions[which(df.summary2$condition == "XY_WT")]
        x_end <- x_positions[which(df.summary2$condition == "XY_Hom")]
        plot <- add_pvalue(plot, x_start, x_end, max(gene_exp$exp) + 0.2 * max(gene_exp$exp), p_XY_Hom)
      }

      # Add p-values for XX comparisons
      if (!is.na(p_XX_Het)) {
        x_start <- x_positions[which(df.summary2$condition == "XX_WT")]
        x_end <- x_positions[which(df.summary2$condition == "XX_Het")]
        plot <- add_pvalue(plot, x_start, x_end, max(gene_exp$exp) + 0.05 * max(gene_exp$exp), p_XX_Het)
      }

      if (!is.na(p_XX_Hom)) {
        x_start <- x_positions[which(df.summary2$condition == "XX_WT")]
        x_end <- x_positions[which(df.summary2$condition == "XX_Hom")]
        plot <- add_pvalue(plot, x_start, x_end, max(gene_exp$exp) + 0.2 * max(gene_exp$exp), p_XX_Hom)
      }

      # Extract legend
      legend <- get_legend(
        plot + theme(legend.direction = "horizontal", legend.box.margin = margin(0, 0, 0, 0))
      )

      # Remove legend from plot
      plot <- plot + theme(legend.position = "none")

      # Store plot
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
        plot.margin = margin(0, 0, 0, 5)
      )
    y_axis <- ggdraw() +
      draw_label(
        "Expression (TPM) vs embryonic stages",
        size = 12,
        hjust = 0.5
      )

    plot_graphs <- plot_grid(plotlist = plotlist, ncol = 3, align = "hv")

    nb_rows_per_type <- ceiling(length(genes)/4)

    plots <- plot_grid(
      gtitle, plot_graphs,
      ncol = 1,
      rel_heights = c(0.1, nb_rows_per_type)
    )
    return(plots)


  } else {

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
          # axis.text.x = element_text(size = 11, angle = 45, hjust=1),
          # axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
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

    plot_graphs <- plot_grid(plotlist = plotlist, ncol = 3, align = "hv")

    nb_rows_per_type <- ceiling(length(genes)/3)

    plots <- plot_grid(
      gtitle, plot_graphs,
      ncol = 1,
      rel_heights = c(0.1, nb_rows_per_type)
    )
    return(plots)
  }
}

##########################################
#                                        #
#          Plot gene expression          #
#                                        #
##########################################

main_nb_rows_per_type <- ceiling(table(markerGenes[markerGenes$Figure == "main", "Cell_type"])/3)
# main_nb_rows_per_type <- main_nb_rows_per_type[unique(markerGenes[markerGenes$Figure == "main", "Cell_type"])]

sup_nb_rows_per_type <- ceiling(table(markerGenes[markerGenes$Figure == "sup", "Cell_type"])/3)
# sup_nb_rows_per_type <- sup_nb_rows_per_type[unique(markerGenes[markerGenes$Figure == "sup", "Cell_type"])]

# main figure
main_plot_genes <- plot_expression_scatter(TPM, markerGenes[markerGenes$Figure == "main", ])

# sup figure
sup_plot_genes <- plot_expression_scatter(TPM, markerGenes[markerGenes$Figure == "sup", ])

# Plot expression of the marker genes
# plot_genes <- plot_expression_scatter(TPM, markerGenes)

main_figure <- plot_grid(
  plotlist = main_plot_genes,
  # labels = "AUTO",
  ncol = 1,
  rel_heights = main_nb_rows_per_type,
  align = "v"
)

sup_figure <- plot_grid(
  plotlist = sup_plot_genes,
  # labels = "AUTO",
  ncol = 1,
  rel_heights = sup_nb_rows_per_type,
  align = "v"
)
##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["main_fig_pdf"]],
  main_figure,
  base_width = 23,
  base_height = 32,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["main_fig_png"]],
  main_figure,
  base_width = 23,
  base_height = 32,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)


# As PDF
save_plot(
  snakemake@output[["sup_fig_pdf"]],
  sup_figure,
  base_width = 23,
  base_height = 32,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["sup_fig_png"]],
  sup_figure,
  base_width = 23,
  base_height = 32,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

# Remove the unwanted file
file.remove("Rplots.pdf")