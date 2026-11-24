# In this document I will retain all of the steps to my most up to date assembly. This will attempt to include relevant decisions made that vary from future steps. 


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#To start I brought in the raw PacBio HiFi data file:

m84082_260311_175030_s3.hifi_reads.bc2024.bam

# Next I had to convert this file to a fastq file. This was the script I used: 

#!/bin/bash

#SBATCH --job-name bam2fastq
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load samtools/1.17

samtools fastq \
m84082_260311_175030_s3.hifi_reads.bc2024.bam \
| gzip > m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

# This sucessfully porducd the fastq file I will use for future input seen below: 

m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#The first step is Genome profiling with Meryl and Genomescope

#Again we start with the raw fastq file below:

m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#Below is the script I ran for Meryl:

#!/bin/bash

#SBATCH --job-name Meryl
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load meryl/1.3

meryl count \
k=31 \
threads=$SLURM_CPUS_PER_TASK \
memory=220 \
output woodrat.k31.meryl \
m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#Next we run Meryl stat to identify the statistics of our Meryl database:

#!/bin/bash

#SBATCH --job-name Merylstat
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load meryl/1.3

meryl statistics woodrat.k31.meryl

#Third, we run Merylhist to generate the histogram for the Meryl database suing the following script:

#!/bin/bash

#SBATCH --job-name Merylhist
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load meryl/1.3

meryl histogram woodrat.k31.meryl > woodrat.k31.hist

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Next we run Genomescope2 in the web browser link below:

http://genomescope.org/genomescope2.0/
  
#I had downloaded my histogram file (woodrat.k31.hist) which I used on the website with the following parameters:
K-mer length
31
Ploidy
2
Max k-mer coverage
-1
Average k-mer coverage for polyploid genome
-1

#This produced plots of the histogram which may be viewed at this link

http://genomescope.org/genomescope2.0/analysis.php?code=MOkGWBikvdnl7V7qww4R

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Next we run Hifiasm in the Hi-C integrated mode. In a new directory  moved over our fastq file and all 4 of our Hi-C data:

m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz
WOODRAT604_L_R1.fastq.gz
WOODRAT604_L_R2.fastq.gz

#To actually run HiFiasm in the integrated mode the script I ran was:

#!/bin/bash

#SBATCH --job-name=HiFiasm
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load hifiasm/0.25.0

hifiasm -o woodrat_hifiandhic.asm -t32 --h1 WOODRAT604_L_R1.fastq.gz --h2 WOODRAT604_L_R2.fastq.gz m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Now we need to QC these runs: First we will run GFAstats: these are the files we start with

woodrat_hifiandhic.asm.hic.hap1.p_ctg.gfa
woodrat_hifiandhic.asm.hic.hap2.p_ctg.gfa

#GFAstats then converst these files into fasta format with this scirpt:

#!/bin/bash

#SBATCH --job-name=gfa
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load gfastats/1.3.6

gfastats \
-f woodrat_hifiandhic.asm.hic.hap1.p_ctg.gfa \
--discover-paths \
-o woodrat_hap1.fasta

gfastats \
-f woodrat_hifiandhic.asm.hic.hap2.p_ctg.gfa \
--discover-paths \
-o woodrat_hap2.fasta

#Next we generate the statistics, first without estimated genome size and then again with that information:

#!/bin/bash

#SBATCH --job-name=stats
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load gfastats/1.3.6

gfastats \
-f woodrat_hap1.fasta \
--stats \
> woodrat_hap1_gfastats.txt

gfastats \
-f woodrat_hap2.fasta \
--stats \
> woodrat_hap2_gfastats.txt

# and with:

#!/bin/bash

#SBATCH --job-name=stats
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load gfastats/1.3.6

gfastats woodrat_hap1.fasta 2742581868 --stats --nstar-report > woodrat_hap1_gfastats_egs.txt

gfastats woodrat_hap2.fasta 2742581868 --stats --nstar-report > woodrat_hap2_gfastats_egs.txt

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Next we ran BUSCO. This started with the following files

woodrat_hap1.fasta
woodrat_hap2.fasta

#This was the shell script:

#!/bin/bash

#SBATCH --job-name=BUSCO
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load busco/5.4.7

busco \
-i woodrat_hap1.fasta \
-m genome \
-l mammalia_odb10 \
-o busco_hap1 \
-c 32

busco \
-i woodrat_hap2.fasta \
-m genome \
-l mammalia_odb10 \
-o busco_hap2 \
-c 32

# To plot busco I took the summary.txt from each haplotype and made a new directory called:

plotdir

# I ran this script to get the BUSCO plots:

#!/bin/bash

#SBATCH --job-name=BUSCO
#SBATCH --partition=shared
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=2:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load busco/5.4.7

generate_plot.py -wd plotdir

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Next we ran Merqury starting with these files and the meryl.db:

woodrat_hap1.fasta
woodrat_hap2.fasta
woodrat.k31.meryl

#This was the shell script:

#!/bin/bash

#SBATCH --job-name=Merqury
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load merqury/1.3

merqury.sh \
woodrat.k31.meryl \
woodrat_hap1.fasta \
woodrat_hap2.fasta \
woodrat_merqury

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Arima mapping pipeline steps: 

#To run the arima mapping pipeline I used the scripts provided by arima modified for my use case. Below I will go through the steps for hap1 first, then repeat for hap2

#First I clones the repository in my directory: 
git clone https://github.com/ArimaGenomics/mapping_pipeline

#Next I rand the bwa index step independently from the rest of the script: 
# STEP 0

#!/bin/bash

#SBATCH --job-name=bwa_index_hap1
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=24GB
#SBATCH --time=5:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw -p bwa_hap1/bwaindex_hap1 woodrat_hap1.fasta

#Step 0:
#Next I am trimming the reads using the following script. This will be done for each of the Hi-C raw files. 
  
#!/bin/bash

#SBATCH --job-name=R1_trim
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=8GB
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17

echo "### Step 0: Trimming the reads" 
zcat WOODRAT_L_R1.fastq.gz | awk '{ if(NR%2==0) {print substr($1,6)} else {print} 
}' | gzip > WOODRAT_L_R1_trimmed.fastq.gz 

#!/bin/bash

#SBATCH --job-name=R1_604_trim
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=8GB
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17

echo "### Step 0: Trimming the reads" 
zcat WOODRAT604_L_R1.fastq.gz | awk '{ if(NR%2==0) {print substr($1,6)} else {print} 
}' | gzip > WOODRAT604_L_R1_trimmed.fastq.gz 

#!/bin/bash

#SBATCH --job-name=R2_trim
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=8GB
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17

echo "### Step 0: Trimming the reads" 
zcat WOODRAT_L_R2.fastq.gz | awk '{ if(NR%2==0) {print substr($1,6)} else {print} 
}' | gzip > WOODRAT_L_R2_trimmed.fastq.gz

#!/bin/bash

#SBATCH --job-name=R2_604_trim
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=8GB
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17

echo "### Step 0: Trimming the reads" 
zcat WOODRAT604_L_R2.fastq.gz | awk '{ if(NR%2==0) {print substr($1,6)} else {print} 
}' | gzip > WOODRAT604_L_R2_trimmed.fastq.gz

#this worked and reduced the raw file sizes by about 15 percent each. NOTE: I will later change the output files to conform with the script

#Now I will run Step 1-3. I will seperate the 3-the end because I need to run 1-3 for each hi-c technical replicate, then merge, then run the last part for each haplotype.
# This is the scirpt I ran for steps 1-3 below. I ensures all raw fastq files were present. Additinally I ran the following script to produce a faidx for the reference:

module load samtools/1.17
samtools faidx woodrat_hap1.fasta

#This produced: 
woodrat_hap1.fasta.fai

#1-3 script:

#!/bin/bash

#SBATCH --job-name=Arima_Mapping
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 11/16/2026 #
##############################################

# Below find the commands used to map HiC data.

# Replace the variables at the top with the correct paths for the locations of files/programs on your system.

# This bash script will map one paired end HiC dataset (read1 & read2 FASTQs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Modules #
##########################################

module load biocontainers/default
module load gcc/11.2.0
module load bwa/0.7.17
module load samtools/1.17
module load picard/3.1.1

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t 12 bwaindex_hap1/bwaindex_hap1 WOODRAT_TRIMMED_R_1.fastq.gz | samtools view -@ 12 -b - > Raw_WOODRAT_TRIMMED_R_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t 12 bwaindex_hap1/bwaindex_hap1 WOODRAT_TRIMMED_R_2.fastq.gz | samtools view -@ 12 -b - > Raw_WOODRAT_TRIMMED_R_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h Raw_WOODRAT_TRIMMED_R_1.bam | perl filter_five_end.pl | samtools view -Sb - > filt_WOODRAT_TRIMMED_R_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h Raw_WOODRAT_TRIMMED_R_2.bam | perl filter_five_end.pl | samtools view -Sb - > filt_WOODRAT_TRIMMED_R_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl two_read_bam_combiner.pl filt_WOODRAT_TRIMMED_R_1.bam filt_WOODRAT_TRIMMED_R_2.bam samtools 10 | samtools view -bS -t woodrat_hap1.fasta.fai - | samtools sort -@ 12 -o Woodrat_hap1.bam -
  
  echo "### Step 3.B: Add read group"
picard AddOrReplaceReadGroups INPUT=Woodrat_hap1.bam OUTPUT=3B_Woodrat_hap1.bam ID=WOODRAT_batch1 LB=WOODRAT SM=WOODRAT PL=ILLUMINA PU=none
###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X # number of the technical replicate set e.g. 1
#	REP_LABEL=${LABEL}_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

This is currently failing at 3A
Test scripts:

#!/bin/bash

#SBATCH --job-name=Arima_Mapping
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 11/16/2026 #
##############################################

# Below find the commands used to map HiC data.

# Replace the variables at the top with the correct paths for the locations of files/programs on your system.

# This bash script will map one paired end HiC dataset (read1 & read2 FASTQs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Modules #
##########################################

module load biocontainers/default
module load gcc/11.2.0
module load bwa/0.7.17
module load samtools/1.17
module load picard/3.1.1

echo "### Step 3A: Pair reads & mapping quality filter"
perl two_read_bam_combiner.pl filt_WOODRAT_TRIMMED_R_1.bam filt_WOODRAT_TRIMMED_R_2.bam ./samtools_wrapper.sh 10 | samtools view -b -t woodrat_hap1.fasta.fai - | samtools sort -@ 12 -o Woodrat_hap1.bam -

echo "### Step 3.B: Add read group"
picard AddOrReplaceReadGroups \
INPUT=Woodrat_hap1.bam \
OUTPUT=3B_Woodrat_hap1.bam \
ID=WOODRAT_batch1 \
LB=WOODRAT \
SM=WOODRAT \ 
PL=ILLUMINA \
PU=none \



  



#From Initial runs in summer
perl two_read_bam_combiner.pl woodrat_hap1_R1_filtered.bam woodrat_hap1_R2_filtered.bam ./samtools_wrapper.sh 10 | samtools view -b -t woodrat_hap1_purged.fa.fai - | samtools sort -@ 128 -o woodrat_hap1_HiC_paired.bam -

# Rest of script:

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X # number of the technical replicate set e.g. 1
#	REP_LABEL=${LABEL}_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 FASTQ file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') # BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"


































