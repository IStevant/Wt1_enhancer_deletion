source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

mapping_samplesheet <- snakemake@input[["samplesheet"]]
count_file <- snakemake@input[["counts"]]
TPM_file <- snakemake@input[["TPM"]]
protein_genes <- snakemake@input[["protein_genes"]]
from_nfcore <- snakemake@params[["from_nfcore"]]
minReads <- snakemake@params[["minReads"]]
minTPM <- snakemake@params[["minTPM"]]
outliers <- snakemake@params[["outliers"]]
remove_outliers <- snakemake@params[["output_files"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Generate the read count matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_gene_matrix <- function(csv_file, genes) {
  # load file
  raw_counts <- read.csv(file = csv_file, header=TRUE, sep = "\t")
  # Remove transcrtpt ID column
  raw_counts <- raw_counts[,-2]
  colnames(raw_counts) <- gsub("[.]", "-", colnames(raw_counts))
  # If matrix come from nf-core, remove unwanted column and remplace gene IDs by gene names
  if(from_nfcore == "yes") {
    # Get gene IDs
    genes <- raw_counts[,1]
    genes <- sapply(strsplit(genes,"\\."), function(x) x[1])
    rownames(raw_counts) <- genes

    # Convert gene ID into gene names
    gene_names <- na.omit(
      as.data.frame(
        AnnotationDbi::mapIds(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = genes,
          column = 'SYMBOL',
          keytype = 'ENSEMBL'
        )
      )
    )
    names(gene_names) <- "gene_symbol"
    merged <- merge(x = gene_names, y = raw_counts, by = 0)
    merged <- merged[,c(-1,-3)]
    merged <- merged[!duplicated(merged$gene_symbol), ]
    rownames(merged) <- merged$gene_symbol
    merged <- merged[,-1]
    colnames(merged) <- gsub("[.]", "-", colnames(merged))

    raw_counts <- merged
  }

  # remove outliers
  if(remove_outliers == "_wo_outliers"){
    raw_counts <- raw_counts[,!colnames(raw_counts) %in% outliers]
  }
  # Transform values as integers for DESeq2 that does not support floats
  raw_counts <- round(raw_counts, digits = 0)

  # Select protein coding genes
  protein_coding_genes <- read.csv(file = protein_genes)
  raw_counts <- raw_counts[rownames(raw_counts) %in% protein_coding_genes$external_gene_name, ]
  raw_counts <- raw_counts[, order(names(raw_counts))]
  return(raw_counts)
}

#' Generate the tpm matrix
#' @param csv_file Path to the read count matrix.
#' @param genes Keep all the genes or only the protein coding genes. Values can be "all" or "protein".
#' @return Return a dataframe.
get_TPM_counts <- function(csv_file, genes) {
  # load file
  tpm <- read.csv(file = csv_file, header=TRUE, sep = "\t")
  # Remove transcrtpt ID column
  tpm <- tpm[,-2]
  colnames(tpm) <- gsub("[.]", "-", colnames(tpm))
  # If matrix come from nf-core, remove unwanted column and remplace gene IDs by gene names
  if(from_nfcore == "yes") {
    # Get gene IDs
    genes <- tpm[,1]
    genes <- sapply(strsplit(genes,"\\."), function(x) x[1])
    rownames(tpm) <- genes

    # Convert gene ID into gene names
    gene_names <- na.omit(
      as.data.frame(
        AnnotationDbi::mapIds(
          org.Mm.eg.db::org.Mm.eg.db,
          keys = genes,
          column = 'SYMBOL',
          keytype = 'ENSEMBL'
        )
      )
    )
    names(gene_names) <- "gene_symbol"
    merged <- merge(x = gene_names, y = tpm, by = 0)
    merged <- merged[,c(-1,-3)]
    merged <- merged[!duplicated(merged$gene_symbol), ]
    rownames(merged) <- merged$gene_symbol
    merged <- merged[,-1]
    colnames(merged) <- gsub("[.]", "-", colnames(merged))

    tpm <- merged
  }

  # remove outliers
  if(remove_outliers == "_clean"){
    tpm <- tpm[,!colnames(tpm) %in% outliers]
  }
  return(tpm)
}

#' Normalize the read count using the size factor normalization from DESeq2
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_normalized_counts <- function(raw_counts, samplesheet) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  # norm_counts <- assay(vst(dds, blind=FALSE))

  return(norm_counts)
}

#' When the maximum expression value (TPM or read count) of a gene between samples is under a certain threshold, we considere it is not properly expressed and the expression values are set to 0.
#' @param row Current gene expression data.
#' @param col_names Save sample names to keep the colnames in the final dataframe.
#' @param minExp Minimum expression value. Default is 10.
#' @return Return a dataframe.
filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}

#' Apply the expression filtration for lowly expressed genes.
#' @param data Read count or TPM matrix.
#' @param minExp Minimum expression value. Default is 10.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp) {
  col_names <- colnames(data)
  data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
  data <- as.data.frame(data)
  data <- data[rowSums(data[]) > 0, ]
  return(data)
}

###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

raw_counts <- get_gene_matrix(
  count_file,
  "protein"
)

TPM <- get_TPM_counts(
  TPM_file,
  "protein"
)


###########################################
#                                         #
#       Filter low expressed genes        #
#                                         #
###########################################

# Remove gene expression if max value < x TPM
TPM <- run_filter_low_counts(TPM, minTPM)
kept_genes <- rownames(TPM[rowSums(TPM) > 0, ])

# Remove gene expression if max value < x reads
raw_counts <- run_filter_low_counts(raw_counts, minReads)

# Remove genes with low TPM from the count matrix
raw_counts <- raw_counts[rownames(raw_counts) %in% kept_genes, ]

###########################################
#                                         #
#             Get samplesheet             #
#                                         #
###########################################

# Generate the samplesheet whith all samples
conditions <- gsub("_Rep.*", "\\1", colnames(raw_counts))
replicate <- sub(".*(Rep[0-9]+)", "\\1", colnames(raw_counts))

samplesheet <- data.frame(
  sample = colnames(raw_counts),
  conditions = conditions,
  replicate = replicate
)

###########################################
#                                         #
#           Normalize read counts         #
#                                         #
###########################################

norm_counts <- get_normalized_counts(
  raw_counts,
  samplesheet
)

###########################################
#                                         #
#             Order matrices              #
#                                         #
###########################################

# If starting from nf-core/RNAseq pipeline, order the matrices according to the sample order in the samplesheet
if (length(mapping_samplesheet)) {
  samplesheet_table <- read.csv(file=mapping_samplesheet, header = TRUE)
  sample_order <- as.vector(samplesheet_table$sample)
  TPM <- TPM[ , match(sample_order, colnames(TPM))]
  raw_counts <- raw_counts[ , match(sample_order, colnames(raw_counts))]
  norm_counts <- norm_counts[ , match(sample_order, colnames(norm_counts))]
}

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

# Save the matrices without ourliers fort the analysis

write.csv(TPM, snakemake@output[["TPM"]])
write.csv(raw_counts, snakemake@output[["counts"]])
write.csv(norm_counts, snakemake@output[["norm_counts"]])
write.csv(samplesheet, snakemake@output[["samplesheet"]], quote = FALSE, row.names = FALSE)
