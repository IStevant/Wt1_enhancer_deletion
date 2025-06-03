source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

cluster_file <- read.csv(snakemake@input[["cluster_file"]], header = TRUE)
sig_DEG_table <-  read.csv(snakemake@input[["sig_DEG_table"]], header = TRUE, sep="\t")

###########################################
#                                         #
#             Merge tables                #
#                                         #
###########################################

# Re,ame columns
colnames(cluster_file) <- c("Genes", "Cluster")

# Merge the tables by the gene names
DEG_merge <- merge(sig_DEG_table, cluster_file, by="Genes")

# Reorder the columns so the cluster ID is just after the padj
DEG_merge <- DEG_merge[, c(1:5, 8, 6,7)]

write.table(DEG_merge, file = snakemake@output[["DEG_table"]], quote = FALSE, row.names = FALSE, sep = "\t")
