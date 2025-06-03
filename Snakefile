'''
Author: Isabelle Stévant
Affiliation: Gonen lab, Bar Ilan university
Date: 02/06/2025
Licence: MIT

The parameters of the analysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the RNA-seq data from nf-code/atacseq pipeline.
'''



# Pipeline configuration file
configfile: "analysis_parameters.yaml"

# Define pipeline output folders
data = f'{config["project_name"]}/data'
output_tables  = f'{config["project_name"]}/tables'
output_png     = f'{config["project_name"]}/graphs/PNG'
output_pdf     = f'{config["project_name"]}/graphs/PDF'

# If starting the analysis from the nf-core/ATACseq output directory
# If not, define the path to the files required to run the downstream analysis
if config["from_nf-core"] == "yes":
    mapping_folder = config["nf-core_path"]
    print("Starting from nf-core/RNAseq output directory: ", mapping_folder)
    samplesheet    = config["nf-core_samplesheet"]
    multiQC        = f'{mapping_folder}/{config["nf-core_multiQC"]}'
    readCounts     = f'{mapping_folder}/{config["nf-core_readCounts"]}'
    TPM            = f'{mapping_folder}/{config["nf-core_TPM"]}'

else:
    print("Starting from custom data directories")
    samplesheet    = []                # create an empty list, will be tested in the R script to know how to order samples
    readCounts     = config["counts"]  # The pre-processed count table
    TPM            = config["TPM"]     # The pre-processed TPM table

# Test if the genome version is defined in the parameter file.
# If not, set "mm10" as default.
try:
    genome = config["genome_version"]
except NameError:
    genome = "mm10"

###########################################
#                                         #
#            Handle outliers              #
#                                         #
###########################################

# If no outliers are declared, generate the expression matrices with all the samples
if config["outliers"] == "None":
    # Ugly trick. 
    # Cannot accept emty string so I removed the "." before extention file to add it with the exp_matrice wildcard...
    exp_matrice = ["."]
else:
# Else if outliers are declared, generate the expression matrices with and without outliers
    exp_matrice = ["_wo_outliers.","_w_outliers."]

###########################################
#                                         #
#             List of output              #
#                                         #
###########################################

# List of output figures
rule_output_list = [
    expand(f"{output_png}/corr_pca{{matrix}}png", matrix=exp_matrice),
    expand(f"{output_png}/DEG_heatmap{{matrix}}png", matrix=exp_matrice),
    expand(f"{output_png}/marker_genes{{matrix}}png", matrix=exp_matrice),
    expand(f"{output_tables}/DEG_with_clusters{{matrix}}tsv", matrix=exp_matrice)
]

# If we start the analysis directly from the nf-core result folder
if config["from_nf-core"] == "yes":
    rule_output_list.append(f"{output_tables}/mapping_summary_count.csv"),
    rule_output_list.append(f"{output_pdf}/RNA_QC_mapping.pdf")

###########################################
#                                         #
#                  Rules                  #
#                                         #
###########################################

# Run the whole pipeline
rule all:
	input: rule_output_list

# Install the necessary R packages using Renv
rule install_packages:
	script:
		"renv/restore.R"

# Get mapping result summary
rule Mapping_summary:
    input: 
        samplesheet = samplesheet,
        multiQC     = multiQC
    params:
        reads       = config["nf-core_reads"]
    output:
        counts      = f"{output_tables}/mapping_summary_count.csv",
        percentages = f"{output_tables}/mapping_summary_percent.csv"
    threads: 1
    resources:
        mem_mb = 4000
    shell: 
        "sh workflow/scripts/RNA_get_mapping_summary.sh \
        {input.samplesheet} \
        {input.multiQC} \
        {params.reads} \
        {output.counts} \
        {output.percentages}"

# Quality control of the data
rule Mapping_report:
    input:
        summary_counts = f"{output_tables}/mapping_summary_count.csv",
        summary_pct    = f"{output_tables}/mapping_summary_percent.csv",
        TPM            = TPM
    params:
        mapping_folder = mapping_folder
    output:
        pdf_report     = f"{output_pdf}/RNA_QC_mapping.pdf",
        png_report     = f"{output_png}/RNA_QC_mapping.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/mapping_qc_plots.R"

# Generate filtered read counts per protein coding genes and TPM matrices
rule Get_matrices:
    input:
        samplesheet   = samplesheet,               # If samplesheet exists, use the order of the samples to order the columns accordingly, if no samplesheet, do nothing
        counts        = readCounts,                # Raw read count per gene matrix comming from the nf-core/rnaseq pipeline
        TPM           = TPM,                       # TPM per gene matrix comming from the nf-core/rnaseq pipeline
        protein_genes = config["protein_genes"]    # List of the mouse protein coding genes
    params:
        from_nfcore   = config["from_nf-core"],    # Precise if the matrices come from nf-core to convert gene IDs to gene names
        minReads      = config["minReads"],        # Minimal number of raw reads from which we consider a gene expressed
        minTPM        = config["minTPM"],          # Minimal number of TPM from which we consider a gene expressed
        outliers      = config["outliers"],        # List of outliers
        output_files  = "{matrix}"                 # Remove or not the outliers
    output:
        TPM           = f"{data}/TPM{{matrix}}csv",            # TPM matrix without the outliers
        counts        = f"{data}/raw_counts{{matrix}}csv",     # Filtered raw read counts
        norm_counts   = f"{data}/norm_counts{{matrix}}csv",    # Filtered normalized read counts (normalization by the library size)
        samplesheet   = f"{data}/samplesheet{{matrix}}csv"     # Description of the samples for downstream analysis
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "workflow/scripts/clean_matrices.R"

# Draw the correlation matrix between samples and the PCA
rule Plot_corr_PCA:
    input:
        norm_counts   = f"{data}/norm_counts{{matrix}}csv",     # Filtered normalized read counts (normalization by the library size)
    params:
        corr_method   = config["corr_met"]                       # Correlation method (can be either "Pearson" or "Spearman")
    output:
        pdf           = f"{output_pdf}/corr_pca{{matrix}}pdf",  # Figure as PDF
        png           = f"{output_png}/corr_pca{{matrix}}png"   # Figure as PNG
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "workflow/scripts/correlation_PCA.R"

# Get genes differentially expressed genes
rule Get_DEGs:
    input:
        raw_counts    = f"{data}/raw_counts{{matrix}}csv",     # Filtered raw read counts
        samplesheet   = f"{data}/samplesheet{{matrix}}csv",    # Description of the samples
        TF_genes      = config["TF_genes"],                              # List of known mouse transcription factors
        TF_pheno      = config["TF_pheno"]                                # List of genes with gonadal phenotypes from MGI OBO database
    params:
        adjpval       = config["adjpval"],    # Adjusted p-value threshold to condider a gene differentially expressed
        log2FC        = config["log2FC"],     # Log2 fold change threshold to condider a gene differentially expressed
    output:
        sig_DEG_obj   = f"{data}/sig_DEGs{{matrix}}Robj",   # R object containing the only the significant DESeq2 results
        sig_DEG_table = f"{data}/DEGs{{matrix}}tsv"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/get_DEG.R"

# Draw the heatmap of the dynamically expressed genes and top GO term enrichment for each cluster side by side.
# Return the genes in each clusters and the GO term enrichment result table.
rule Plot_heatmap_GO:
    input:
        TF_genes     = config["TF_genes"],                              # List of known mouse transcription factors
        sig_DEGs     = f"{data}/sig_DEGs{{matrix}}Robj",     # Robj with the filtered dynamic genes
        norm_counts  = f"{data}/norm_counts{{matrix}}csv",    # Filtered normalized read counts (normalization by the library size)
        samplesheet  = f"{data}/samplesheet{{matrix}}csv",    # Description of the samples
        marker_genes = config["marker_genes"],                          # List of known mouse transcription factors
        TF_pheno     = config["TF_pheno"]                               # List of genes with gonadal phenotypes from MGI OBO database
    output:
        GO           = f"{output_tables}/GO_DEG{{matrix}}csv",                  # Simplified GO term enrichment result table
        cluster_file = f"{data}/DEG_heatmap_clusters{{matrix}}csv",    # Genes per clusters
        pdf          = f"{output_pdf}/DEG_heatmap{{matrix}}pdf",
        png          = f"{output_png}/DEG_heatmap{{matrix}}png"
    threads: 1
    resources:
        mem_mb = 16000
    script:
        "workflow/scripts/DEG_heatmap_GO.R"

# Add the cluster id to the DE gene result table
rule DEG_result_table:
    input:
        cluster_file = f"{data}/DEG_heatmap_clusters{{matrix}}csv",
        sig_DEG_table = f"{data}/DEGs{{matrix}}tsv"
    output:
        DEG_table = f"{output_tables}/DEG_with_clusters{{matrix}}tsv"
    threads: 1
    resources:
        mem_mb = 16000
    script:
        "workflow/scripts/DEG_result_table_with_clusters.R"


# Get genes differentially expressed genes
rule Plot_marker_genes:
    input:
        TPM          = f"{data}/TPM{{matrix}}csv",
        marker_genes = config["marker_genes"]
    output:
        fig_png = f"{output_png}/marker_genes{{matrix}}png",
        fig_pdf = f"{output_pdf}/marker_genes{{matrix}}pdf",
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/plot_marker_genes.R"

onsuccess:
    from snakemake.report import auto_report
    auto_report(workflow.persistence.dag, f'{config["project_name"]}/report/report.html')
    