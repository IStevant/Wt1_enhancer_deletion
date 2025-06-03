source(".Rprofile")

###########################################
#                                         #
#             Load libraries              #
#                                         #
###########################################

library("ggplot2")
library("scales")
library("cowplot")
library("shades")
library("dplyr")

###########################################
#                                         #
#                Load data                #
#                                         #
###########################################

summary_counts <- snakemake@input[["summary_counts"]]
summary_pct <- snakemake@input[["summary_pct"]]
TPM <- snakemake@input[["TPM"]]

samples <- as.vector(read.table(summary_counts, header=TRUE, sep=",")[,1])


# Get the summary tables for the samples specified in the sample list 
get_table <- function(sample){

	summary_table_count <- read.table(summary_counts, header=TRUE, sep=",")
	summary_table_count <- summary_table_count[grep(sample, summary_table_count$sample),]

	summary_table_pct <- read.table(summary_pct, header=TRUE, sep=",")
	summary_table_pct <- summary_table_pct[grep(sample, summary_table_pct$sample),]

	return(list(counts=summary_table_count, pct=summary_table_pct))
}

summaries <- lapply(samples, function(sample) get_table(sample))

count_summary <- data.table::rbindlist(lapply(summaries, function(tables) tables$counts))
pct_summary <- data.table::rbindlist(lapply(summaries, function(tables) tables$pct))
colnames(pct_summary) <- stringr::str_replace(colnames(pct_summary), "X..", "%.")

mapping_summary <- cbind(count_summary, pct_summary[,!"sample"])
colnames(mapping_summary) <-   stringr::str_to_sentence(colnames(mapping_summary))

samples <- mapping_summary$Sample

conditions <- gsub("_Rep.*", "\\1", samples)

group <- unique(conditions)

mapping_summary <- data.frame(
	group=conditions,
	mapping_summary
)

# Pick colors
set.seed(1234)
colours <- randomcoloR::distinctColorPalette(length(group))
names(colours) <- group

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################


plot_reads <- function(mapping_summary){
	print("Plot total reads")

	data <- mapping_summary
	data$label <- paste0(round(data$Total/1000000, 1), "M")
	data$Total <- data$Total/1000000
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	reads <- ggplot(data, aes(Total, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Million read counts") +
		theme_light() +
		ggtitle("Read counts per samples") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(reads)
}

plot_mapped <- function(mapping_summary){
	print("Plot mapped reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Mapped, "%")
	data$Mapped <- data$Mapped/1000000
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	mapping <- ggplot(data, aes(X..Mapped, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentage") +
		theme_light() +
		ggtitle("% of mapped reads per samples") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(mapping)
}

plot_duplicated <- function(mapping_summary){
	print("Plot duplicated reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Duplicated, "%")
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	dup <- ggplot(data, aes(X..Duplicated, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of duplicated reads (rel. to mapped reads)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

plot_multimapped <- function(mapping_summary){
	print("Plot multimapped reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Multimapped, "%")
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	dup <- ggplot(data, aes(X..Multimapped, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of multimapped reads (rel. to mapped reads)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

plot_exonic <- function(mapping_summary){
	print("Plot exonic reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Exonic, "%")
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	dup <- ggplot(data, aes(X..Exonic, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of exonic reads (rel. to mapped reads)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

plot_intronic <- function(mapping_summary){
	print("Plot intronic reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Intronic, "%")
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	dup <- ggplot(data, aes(X..Intronic, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of intronic reads (rel. to mapped reads)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

plot_intergenic <- function(mapping_summary){
	print("Plot intergenic reads")

	data <- mapping_summary
	data$label <- paste0(data$X..Intergenic, "%")
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	dup <- ggplot(data, aes(X..Intergenic, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of intergenic reads (rel. to mapped reads)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

plot_nb_genes <- function(tpm){
	print("Plot Nb of detected genes")

	TPM <- read.csv(tpm, sep="\t", row.names=1)
	TPM <- TPM[,-1]
	colnames(TPM) <- gsub("[.]", "-", colnames(TPM))

	TPM[TPM>0] <- 1

	gene_per_sample <- colSums(TPM)

	data <- data.frame(
		Sample = names(gene_per_sample),
		Gene_count = as.vector(gene_per_sample),
		group=gsub("_Rep.*", "\\1", colnames(TPM))
	)

	data$label <- scales::comma(data$Gene_count)
	ratio <- 1
	fontsize <- 8
	text_size <- 2.5
	legend_position <- "none"

	print(names(gene_per_sample))
	print(mapping_summary$Sample)

	dup <- ggplot(data, aes(Gene_count, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="#000000", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Gene count") +
		theme_light() +
		ggtitle("Nb of detected genes") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(dup)
}

reads <- plot_reads(mapping_summary)
mapping <- plot_mapped(mapping_summary)
duplicate <- plot_duplicated(mapping_summary)
multimapped <- plot_multimapped(mapping_summary)
exonic <- plot_exonic(mapping_summary)
intronic <- plot_intronic(mapping_summary)
intergenic <- plot_intergenic(mapping_summary)
genes <- plot_nb_genes(TPM)

figure <- plot_grid(
	reads,
	mapping,
	duplicate,
	multimapped,
	exonic,
	intronic,
	intergenic,
	genes,
	ncol=2, 
	align="hv", 
	axis="bt", 
	rel_widths = c(1, 1),
	labels="AUTO",
	label_size = 8
)

save_plot(
	snakemake@output[["pdf_report"]], 
	figure, 
	dpi=300, 
	bg = "white", 
	base_width=20, 
	base_height=32,
	units = c("cm")
)

save_plot(
	snakemake@output[["png_report"]], 
	figure, 
	dpi=300, 
	bg = "white", 
	base_width=20, 
	base_height=35,
	units = c("cm")
)

# Remove the unwanted file
file.remove("Rplots.pdf")