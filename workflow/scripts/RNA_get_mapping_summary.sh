#!/bin/bash


# Usage: sh mapping_summary.sh samplesheet.csv mapping_result_folder single|paired output_counts output_percent

# Get sample name and final result file name
samples=(`awk -F "\"*,\"*" 'NR>1{print $1}' $1`)
results=$2
res_count_name=$4
res_pct_name=$5

# If the RNA-seq is single-end
if [[ "$3" == *"single"* ]] ; then
	# Column names for count table
	echo -e "sample,total,mapped,duplicated,multimapped,exonic,intronic,intergenic" > $res_count_name
	# Column names for percentage table
	echo -e "sample,% mapped,% duplicated,% multimapped,% exonic,% intronic,% intergenic" > $res_pct_name

	# For each sample
	for sample in "${samples[@]}"; do
		# Get the total number of reads in fastq files
		tot_reads=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $4}' | sed 's/\.0//g'`
		echo $tot_reads
		# Mapped reads
		mapped=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $8}' | sed 's/\.0//g'`
		pct_mapped=`awk -v a=$mapped -v b=$tot_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Duplicated reads
		dup=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $13}' | sed 's/\.0//g'`
		pct_dup=`awk -v a=$dup -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Multi mapped reads
		multi=`grep $sample $results/multiqc_rsem.txt | awk '{print $8}' | sed 's/\.0//g'`
		pct_multi=`awk -v a=$multi -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Exonic reads
		exonic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $2}' | sed 's/\.0//g'`
		pct_exonic=`awk -v a=$exonic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Intronic reads
		intronic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $3}' | sed 's/\.0//g'`
		pct_intronic=`awk -v a=$intronic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Intergenic reads
		intergenic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $4}' | sed 's/\.0//g'`
		pct_intergenic=`awk -v a=$intergenic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		
		sum_feature_reads=$(( exonic+intronic+intergenic ))
		pct_exonic=`awk -v a=$exonic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		pct_intronic=`awk -v a=$intronic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		pct_intergenic=`awk -v a=$intergenic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`


		counts=`echo -e "$sample,$tot_reads,$mapped,$dup,$multi,$exonic,$intronic,$intergenic"`
		percentages=`echo -e "$sample,$pct_mapped,$pct_dup,$pct_multi,$pct_exonic,$pct_intronic,$pct_intergenic"`

		echo $counts >> $res_count_name
		echo $percentages >> $res_pct_name
	done
else
		# Column names for count table
	echo -e "sample,total (R1+R2),paired (R1+R2),mapped (R1+R2),duplicated (R1+R2),multimapped (R1+R2),exonic (R1+R2),intronic (R1+R2),intergenic (R1+R2)" > $res_count_name
	# Column names for percentage table
	echo -e "sample,% paired,% mapped,% duplicated,% multimapped,% exonic,% intronic,% intergenic" > $res_pct_name
	for sample in "${samples[@]}"; do
		# Get the total number of reads in fastq files
		tot_reads=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $4}' | sed 's/\.0//g'`
		echo $tot_reads
		# Get the total number of reads in fastq files
		tot_paired_reads=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $12}' | sed 's/\.0//g'`
		pct_paired=`awk -v a=$tot_paired_reads -v b=$tot_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		echo $tot_paired_read
		# Mapped reads
		mapped=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $8}' | sed 's/\.0//g'`
		pct_mapped=`awk -v a=$mapped -v b=$tot_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Duplicated reads
		dup=`grep $sample $results/multiqc_samtools_stats.txt | awk '{print $13}' | sed 's/\.0//g'`
		pct_dup=`awk -v a=$dup -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Multi mapped reads
		multi=`grep $sample $results/multiqc_rsem.txt | awk '{sum+=$8;} END{print sum;}' | sed 's/\.0//g'`
		pct_multi=`awk -v a=$multi -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Exonic reads
		exonic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $2}' | sed 's/\.0//g'`
		pct_exonic=`awk -v a=$exonic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Intronic reads
		intronic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $3}' | sed 's/\.0//g'`
		pct_intronic=`awk -v a=$intronic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		# Intergenic reads
		intergenic=`grep $sample $results/mqc_qualimap_genomic_origin_1.txt | awk '{print $4}' | sed 's/\.0//g'`
		pct_intergenic=`awk -v a=$intergenic -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		
		sum_feature_reads=$(( exonic+intronic+intergenic ))
		pct_exonic=`awk -v a=$exonic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		pct_intronic=`awk -v a=$intronic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`
		pct_intergenic=`awk -v a=$intergenic -v b=$sum_feature_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`


		counts=`echo -e "$sample,$tot_reads,$tot_paired_read,$mapped,$dup,$multi,$exonic,$intronic,$intergenic"`
		percentages=`echo -e "$sample,$pct_paired,$pct_mapped,$pct_dup,$pct_multi,$pct_exonic,$pct_intronic,$pct_intergenic"`

		echo $counts >> $res_count_name
		echo $percentages >> $res_pct_name

	done
fi
