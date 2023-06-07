#!/bin/bash

# PIPELINE FOR VARIANT CALLING USING GATK4 

# by Eliana Buenaventura (elianabuenaventura@gmail.com)



# PRE/PROCESSING OF RAW READS ############################################################


## Prepare your folder structure #####################################

## Here is the folder structure inside the folder genepanels_gatk. Please note that each
## run has its own folder with the exact same folder structure inside. Inside the folder
## genepanels_gatk, we should also have a folder called known_sites with our truth sets.

   
├── genepanels_gatk
│   ├── known_sites
│   ├── run9
│   │   ├── data
│   │   │   ├── bam
│   │   │   ├── untrimmed_fastq
│   │   │   └── vcf
│   │   ├── docs
│   │   ├── report
│   │   └── results
│   │       ├── coverage
│   │       └── metadata
│   └── run10
│       ├── data
│       │   ├── bam
│       │   ├── untrimmed_fastq
│       │   └── vcf
│       ├── docs
│       ├── report
│       └── results
│           ├── coverage
│           └── metadata
├── ref_genome
├── runs
├── software
└── targets_bed


## Most of the scripts in this pipeline should be run from your particular run folder. Exmaple, from /run9 or /run10.



## Generate an unmapped BAM from FASTQ ################################

## Step 1: Convert FASTQ to uBAM and add read group information using FastqToSam ###
## Make sure your READ_GROUP_NAME, LIBRARY_NAME and PLATFORM_UNIT are adjusted to your data.

cd .../genepanels_gatk/run9

for infile in data/untrimmed_fastq/*_R1_001.fastq.gz
do
	base=$(basename ${infile} _R1_001.fastq.gz)
   java -Xmx8G -jar /mnt/d2/genepanels/software/picard/picard.jar FastqToSam \
	FASTQ=${infile} \
    FASTQ2=data/untrimmed_fastq/${base}_R2_001.fastq.gz \
    OUTPUT=data/bam/${base}_fastqtosam.bam \
    READ_GROUP_NAME=000000000-KRNF3.1 \
    SAMPLE_NAME=${base} \
    LIBRARY_NAME=20230224-123656 \
    PLATFORM_UNIT=000000000-KRNF3.1.20230224-123656 \
    PLATFORM=illumina \
    SEQUENCING_CENTER=AHH-KBA
    echo ${base}
done



## Step 2: Mark adapter sequences using MarkIlluminaAdapters ###

for infile in data/bam/*_fastqtosam.bam
do
	base=$(basename ${infile} _fastqtosam.bam)
   java -Xmx8G -jar /mnt/d2/genepanels/software/picard/picard.jar MarkIlluminaAdapters \
   I=${infile} \
   O=data/bam/${base}_markilluminaadapters.bam \
   M=data/bam/${base}_markilluminaadapters_metrics.txt
done



## Step 3: Align reads with BWA-MEM and merge with uBAM using MergeBamAlignment ###

for infile in data/bam/*_markilluminaadapters.bam
do
	base=$(basename ${infile} _markilluminaadapters.bam)
   java -Xmx8G -jar /mnt/d2/genepanels/software/picard/picard.jar SamToFastq \
   I=${infile} \
   FASTQ=/dev/stdout \
   CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \

   bwa mem -M -t 7 -p /mnt/d2/genepanels/ref_genome/hg38-v0-Homo_sapiens_assembly38.fasta \
   /dev/stdin | \

   java -Xmx16G -jar /mnt/d2/genepanels/software/picard/picard.jar MergeBamAlignment \
   R=/mnt/d2/genepanels/ref_genome/hg38-v0-Homo_sapiens_assembly38.fasta \
   UNMAPPED_BAM=data/bam/${base}_fastqtosam.bam \
   ALIGNED_BAM=/dev/stdin \
   O=data/bam/${base}_piped.bam \
   CREATE_INDEX=true \
   ADD_MATE_CIGAR=true \
   CLIP_ADAPTERS=false \
   CLIP_OVERLAPPING_READS=true \
   INCLUDE_SECONDARY_ALIGNMENTS=true \
   MAX_INSERTIONS_OR_DELETIONS=-1 \
   PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
   ATTRIBUTES_TO_RETAIN=XS
done



# MARK DUPLICATES ########################################################################

for infile in data/bam/*_piped.bam
do
	base=$(basename ${infile} _piped.bam)
   java -Xmx8G -jar /mnt/d2/genepanels/software/picard/picard.jar MarkDuplicates \
   I=${infile} \
   O=data/bam/${base}_marked_duplicates.bam \
   M=data/bam/${base}_marked_dup_metrics.txt
done



# BASE (QUALITY SCORE) RECALIBRATION #####################################################

for infile in data/bam/*_marked_duplicates.bam
do
	base=$(basename ${infile} _marked_duplicates.bam)
   gatk --java-options -Xmx8G BaseRecalibrator \
   -I ${infile} \
   -R /mnt/d2/genepanels/ref_genome/hg38-v0-Homo_sapiens_assembly38.fasta \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-1000G_omni2.5.hg38.vcf.gz \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-hapmap_3.3.hg38.vcf.gz \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites /mnt/d2/genepanels/genepanels_gatk/known_sites/hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -O data/bam/${base}_recal_data.table
done



for infile in data/bam/*_marked_duplicates.bam
do
	base=$(basename ${infile} _marked_duplicates.bam)
   gatk --java-options -Xmx8G ApplyBQSR \
   -I ${infile} \
   -R /mnt/d2/genepanels/ref_genome/hg38-v0-Homo_sapiens_assembly38.fasta \
   --bqsr-recal-file data/bam/${base}_recal_data.table \
   -O data/bam/${base}_BQSR_marked_duplicates.bam
done



# CALL VARIANTS PER-SAMPLE ###############################################################

for infile in data/bam/*_BQSR_marked_duplicates.bam
do
	base=$(basename ${infile} _BQSR_marked_duplicates.bam)
   gatk --java-options -Xmx8G HaplotypeCaller \
   -R /mnt/d2/genepanels/ref_genome/hg38-v0-Homo_sapiens_assembly38.fasta \
   -I ${infile} \
   -O data/bam/${base}_BQSR_marked_duplicates.vcf
done



# VARIANT FILTERING BY HARD FILTERING ###################################################


## Step 1: Select SNPs and indels


### Code - Selects SNPs for several samples

for infile in data/bam/*_BQSR_marked_duplicates.vcf
do
	base=$(basename ${infile} .vcf)
   gatk --java-options -Xmx8G SelectVariants \
      -select-type SNP \
      -V ${infile} \
      -O data/bam/${base}_snps.vcf
done



### Code - Selects indels for several samples

for infile in data/bam/*_BQSR_marked_duplicates.vcf
do
	base=$(basename ${infile} .vcf)
   gatk --java-options -Xmx8G SelectVariants \
      -select-type INDEL \
      -V ${infile} \
      -O data/bam/${base}_indels.vcf
done



## Step 2: Filtering of variants with VariantFiltration.


### Code - Filters SNPs on several samples in same folder - QD < 4.0

for infile in data/bam/*_snps.vcf
do
	base=$(basename ${infile} _snps.vcf)
   gatk --java-options -Xmx8G VariantFiltration \
    -V ${infile} \
    -filter "QD < 4.0" --filter-name "QD4" \
    -filter "QUAL < 32.0" --filter-name "QUAL32" \
    -filter "SOR > 6.0" --filter-name "SOR6" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O data/bam/${base}_snps_filtered.vcf
done



### Code - Filters indels on several samples in same folder - QD < 7.5

for infile in data/bam/*_indels.vcf
do
	base=$(basename ${infile} _indels.vcf)
   gatk --java-options -Xmx8G VariantFiltration \
    -V ${infile} \
    -filter "QD < 7.5" --filter-name "QD7.5" \
    -filter "QUAL < 80.0" --filter-name "QUAL80" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O data/bam/${base}_indels_filtered.vcf
done



## Step 3: Concatenate SNPs and indels in one VCF file

for infile in data/bam/*_snps_filtered.vcf
do
   base=$(basename ${infile} _BQSR_marked_duplicates_snps_filtered.vcf)
   vcf-concat ${infile} data/bam/${base}_BQSR_marked_duplicates_indels_filtered.vcf > data/bam/${base}_var_filtered.vcf
   echo ${base}
done



## Step 4: Select variants in specific regions. 


### Code - selects FMF-pakke targets for several samples (tested)

for infile in data/bam/mefv*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   vcftools --vcf ${infile} --bed /mnt/d2/genepanels/targets_bed/FMF_pakke_targets_hg38.bed \
   --recode --recode-INFO-all -c > data/vcf/${base}_fmf_var_filtered.vcf
   echo ${base}
done



### Code - selects FHH-ADH-pakke targets for several samples (tested)

for infile in data/bam/fhh*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   vcftools --vcf ${infile} --bed /mnt/d2/genepanels/targets_bed/FHH-ADH_pakke_targets_hg38.bed \
   --recode --recode-INFO-all -c > data/vcf/${base}_fhh-adh_var_filtered.vcf
   echo ${base}
done



# CREATE FINAL VCF FILE WITH ONLY PASSING VARIANTS #####################################


cd /mnt/d2/genepanels/genepanels_gatk/runX

for infile in data/vcf/*.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   vcftools --vcf ${infile} --remove-filtered-all \
   --recode --recode-INFO-all -c > data/vcf/${base}_final.vcf
   echo ${base}
done


# VARIANT REPORTING #####################################################################


## PREPARE VARIANT METADATA FOR REPORT WITH VEP #########################################


## Step 1: Produce all needed metadata for variants in coding and non-coding regions for 
## internal reference. Output in table format.

cd .../ensembl-vep

for infile in .../genepanels_gatk/run9/data/vcf/*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   ./vep -i ${infile} --cache \
   --force_overwrite --check_existing --hgvs --sift b --canonical --symbol --numbers --tab \
   --stats_file /mnt/d2/genepanels/genepanels_gatk/run14v/results/metadata/${base}_statsfile.txt \
   --fields "Uploaded_variation,Location,Allele,SYMBOL,Feature,Feature_type,CANONICAL,Consequence,INTRON,EXON,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SIFT,HGVSc,HGVSp,CLIN_SIG" \
   -o STDOUT | grep -v "##" > /mnt/d2/genepanels/genepanels_gatk/run14v/results/metadata/${base}_metadata.txt
   echo ${base}
done



## Step 2: Produce only needed metadata for variants in coding and non-coding regions for 
## internal report. Output in table format.

cd .../ensembl-vep

for infile in .../genepanels_gatk/run9/data/vcf/*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   ./vep -i ${infile} --cache \
   --force_overwrite --check_existing --hgvs --sift b --canonical --symbol --numbers --tab \
   --stats_file /mnt/d2/genepanels/genepanels_gatk/run14/results/metadata/${base}_statsfile_report.txt \
   --fields "Uploaded_variation,Allele,CANONICAL,Consequence,INTRON,EXON,Amino_acids,Codons,Existing_variation,HGVSc,HGVSp,CLIN_SIG" \
   -o STDOUT | ./filter_vep --filter "CANONICAL is YES" | grep -v "##" > /mnt/d2/genepanels/genepanels_gatk/run14/results/metadata/${base}_metadata_report.txt
   echo ${base}
done


## Step 3: Prepare metadata table for internal report. Remove # from table and edit headers.

for infile in results/metadata/*_metadata_report.txt
do
    base=$(basename ${infile} _metadata_report.txt)
    awk -F"\t" 'BEGIN{OFS="\t"} \
    {if(NR==1){$1="Variant_ID"} \
    if(NR==1){$5="Intron"} \
    if(NR==1){$6="Exon"} \
    if(NR==1){$7="AA"} \
    if(NR==1){$12="Clin_Sig"} \
    ; print $1, $2, $4, $5, $6, $7}' \
    $infile > report/${base}_var_metadata_report_p1.txt
done



for infile in results/metadata/*_metadata_report.txt
do
    base=$(basename ${infile} _metadata_report.txt)
    awk -F"\t" 'BEGIN{OFS="\t"} \
    {if(NR==1){$1="Variant_ID"} \
    if(NR==1){$5="Intron"} \
    if(NR==1){$6="Exon"} \
    if(NR==1){$7="AA"} \
    if(NR==1){$12="Clin_Sig"} \
    ; print $1, $9, $10, $11, $12}' \
    $infile > report/${base}_var_metadata_report_p2.txt
done



## PREPARE VARIANTS FOR REPORT ##########################################


## Step 1: Variants from VCF format to table format not for report but only for 
## internal reference.

for infile in data/vcf/*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   gatk VariantsToTable \
     -V ${infile} \
     -F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL \
     -F DP -F QD -F FS -F SOR -F MQRankSum -F ReadPosRankSum \
     -GF AD -GF DP -GF GT -GF GQ -GF PL \
     -O report/${base}_var_filtered_fulltable.table
done



## Step 2: Variants from VCF format to table format for internal report.

for infile in data/vcf/*_var_filtered.vcf
do
   base=$(basename ${infile} _var_filtered.vcf)
   gatk VariantsToTable \
     -V ${infile} \
     -F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL \
     -F DP -F QD -F SOR -GF GT -GF DP -GF GQ -GF AD \
     -O report/${base}_var_filtered.table
done



## Step 3: Table formatting not for report but only for internal reference.

for infile in report/*_var_filtered_fulltable.table
do
    base=$(basename ${infile} _var_filtered_fulltable.table)
    awk -F"\t" 'BEGIN{OFS="\t"} \
    {if(NR==1){$6="Quality (good if > 30)"} \
    if(NR==1){$7="Read Depth (good if > 50)"} \
    if(NR==1){$8="Quality by Depth (good if > 2.0)"} \
    if(NR==1){$9="Fisher Strand (good if < 60)"} \
    if(NR==1){$10="Strand Odd Ratio (good if < 3.0)"} \
    if(NR==1){$11="MQRankSum (good if > -12.5)"} \
    if(NR==1){$12="ReadPosRankSum (good if > -8.0)"} \
    if(NR==1){$13="Allelic Depths"} \
    if(NR==1){$14="Genotype Read Depth"} \
    if(NR==1){$15="Genotype"} \
    if(NR==1){$16="Genotype Quality"} \
    if(NR==1){$17="Phred-scaled Likelihood"} \
    ; print $0}' \
    $infile > report/${base}_var_filtered_fulltable_read.table
done



## Step 4: Table formatting for internal report.

for infile in report/*_var_filtered.table
do
    base=$(basename ${infile} _filtered.table)
    awk -F"\t" 'BEGIN{OFS="\t"} \
    {if(NR==1){$7="RD"} \
    if(NR==1){$8="QD"} \
    if(NR==1){$9="SOR"} \
    if(NR==1){$10="GT"} \
    if(NR==1){$11="GT_RD"} \
    if(NR==1){$12="GT_Q"} \
    if(NR==1){$13="AD"} \
    ; print $0}' \
    $infile > report/${base}_var_filtered_pdf.table
done



## PREPARE COVERAGE INFO FOR REPORT ##########################################

## Step 1: Obtain read depth

### Code - for several samples and Fever-pakke targets

for infile in data/bam/mefv*_BQSR_marked_duplicates.bam
do
	base=$(basename ${infile} _BQSR_marked_duplicates.bam)
	samtools depth -d 0 -a ${infile} -b /mnt/d2/genepanels/targets_bed/Fever_pakke_targets_hg38.bed \
	-o results/coverage/${base}_fever_gatk_dep.cov
	echo ${base}
done



### Code - for several samples and FHH-ADH-pakke targets

for infile in data/bam/casr*_BQSR_marked_duplicates.bam
do
	base=$(basename ${infile} _BQSR_marked_duplicates.bam)
	samtools depth -d 0 -a ${infile} -b /mnt/d2/genepanels/targets_bed/FHH-ADH_pakke_targets_hg38.bed \
	-o results/coverage/${base}_fhh-adh_gatk_dep.cov
	echo ${base}
done



## Step 2: Calculate quality of read depth 

### Code for several samples - threshold 50

for infile in results/coverage/*_gatk_dep.cov
do
	base=$(basename ${infile} _gatk_dep.cov)
	awk 'BEGIN{chr="nothing"}{if($3<50){if($1!=chr||$2!=(end+1)){if(chr!="nothing"&&flag==1){print chr"\t"start"\t"end"\t"cov/len"\t"len ;flag=0};};if(flag==0){chr=$1;start=$2;end=$2;cov=$3;len=1;}else{end=$2;cov=cov+$3;len=len+1;};flag=1;};}END{if(flag==1){print chr"\t"start"\t"end"\t"cov/len"\t"len}}' ${infile} > results/coverage/${base}_gatk_Q50dep.txt
	echo ${base}
done



## Step 3: Format output file

cd results/coverage/

for infile in *_gatk_Q50dep.txt
do
    base=$(basename ${infile} _gatk_Q50dep.txt)
    awk "{print \"${base}\t 	 \" \$0}" "$infile" > /dev/stdout | \
    awk 'BEGIN{print "Sample\tChr\tStart\tEnd\tMean_read_depth\tLength"}1' /dev/stdin >> ../../report/${base}_gatk_belowQ50dep.txt
    echo ${base}
done



# END ####################################################################################
