source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

count <- read.csv(file = snakemake@input[["raw_counts"]], row.names = 1)
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], header = TRUE)
TF_list <- read.csv(file = snakemake@input[["TF_genes"]], header = FALSE)$V1
pheno_TFs <- read.table(file = snakemake@input[["TF_pheno"]], header = FALSE, sep = "\t")
adj.pval <- snakemake@params[["adjpval"]]
log2FC <- snakemake@params[["log2FC"]]

###########################################
#                                         #
#            DESeq analysis               #
#                                         #
###########################################


dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count,
  colData = samplesheet,
  design = ~conditions
)

ddsTC <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
resTC <- DESeq2::results(ddsTC)
resTC$symbol <- GenomicRanges::mcols(ddsTC)$symbol


filtered_DEGs <- subset(resTC, padj < adj.pval)
filtered_DEGs <- subset(filtered_DEGs, abs(log2FoldChange) > log2FC)

filtered_DEGs <- filtered_DEGs[, !colnames(filtered_DEGs) %in% c("lfcSE", "stat")]


filtered <- data.frame(
  filtered_DEGs,
  is.TF = ifelse(rownames(filtered_DEGs) %in% TF_list, "Yes", "-")
)

rownames(pheno_TFs) <- pheno_TFs[, 1]
colnames(pheno_TFs) <- c("genes", "Phenotype")
merged_filtered <- merge(filtered, pheno_TFs[, 2, drop = FALSE], by = 0, all.x = TRUE)
merged_filtered[is.na(merged_filtered)] <- "-"
names(merged_filtered)[names(merged_filtered) == "Row.names"] <- "Genes"
rownames(merged_filtered) <- merged_filtered$Genes

filtered_DEGs <- merged_filtered

write.table(filtered_DEGs[order(filtered_DEGs$padj), ], file = snakemake@output[["sig_DEG_table"]], quote = FALSE, row.names = FALSE, sep = "\t")

filtered_DEGs <- rownames(filtered_DEGs)

save(filtered_DEGs, file = snakemake@output[["sig_DEG_obj"]])
