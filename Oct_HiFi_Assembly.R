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

#The first step is Genome profile 