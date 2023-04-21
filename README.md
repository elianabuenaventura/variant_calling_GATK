# variant_calling_GATK

## pipeline_vc_gatk.sh is a collection of scripts for variant calling from NGS sequence data using GATK 
This is an analysis pipeline (based on GATK, Picard, samtools, bedtools, bcftools, vcftools, BWA, ensembl-vep) for analyzing data collected from a gene panel on a local supercomputer.

This is a wrapper shell script for automating the pipeline starting from pre-processing reads and ending with reporting variants. 
Pipeline runs in loops designed for processing batches of samples. 

It includes scripts that allow for:
* read pre-processing with Picard
* read alignment with BWA and Picard (using GRCh38/hg38 as reference)
* variant calling with HaplotypeCaller of GATK
* annotating and reporting relevant variants with ensembl-vep, GATK and IGV

This version has been tested only on a local supercomputer.
