#Below I have documented my steps assembling the woodrats genome. Notes will include every iteration of each step.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Step 1: PacBio HiFi Data: Adapter filtering using HiFiAdapterFilt
#Documentation: https://github.com/sheinasim-USDA/HiFiAdapterFilt

# To run HiFiAdapterFilt in anvil I first cloned the github directory to my WD using the code below. Upon completion I had a new directory named HiFiAdapterFilt which is called
# In the data filtering step
git clone https://github.com/sheinasim/HiFiAdapterFilt.git

#To ensure HiFiAdapterFilt runs properly I needed to download the conda environment from the ccgp github which has the ncbi database and dependencies to run blast/2.6
#To download this conda environment I ran:
module load anaconda/2024.02-py311

wget https://raw.githubusercontent.com/ccgproject/ccgp_assembly/main/workflows/conda_env/conda.env.ncbi26.yml #This downloads the conda environment directly to my WD

conda env create -f conda.env.ncbi26.yml #This actually creates the conda environment to our WD. 
#This entire step was later rendered usless.

#With all the dependencies loaded in this conda we can now run the shell scirpt below to actually run HiFiAdapterFilt


#!/bin/bash

#SBATCH --job-name hififilt
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=4:00:00

cd $SLURM_SUBMIT_DIR

module load anaconda/2024.02-py311

source $(conda info --base)/etc/profile.d/conda.sh
conda activate ncbi26

TOOL_DIR=/anvil/scratch/x-srudnick/Woodrat_WGS/HiFiAdapterFilt

bash "$TOOL_DIR/hifiadapterfilt.sh" -p m84082_260311_175030_s3.hifi_reads.bc2024 -t 8

#I had issues getting the conda to run. Initially I waswriting the script this was so that the code was the exact same as the CCGP pipeline inducing Blast+. I later identified that
#Blast+ is available in anvil under the name Blast-plus. This allowed me to pivot from the conda environment and simply utilize HiFiAdapterFilt with the most updated packages
#Below is the true steps I used to produce HiFi Filtering:

# To run HiFiAdapterFilt in anvil I first cloned the github directory to my WD using the code below. Upon completion I had a new directory named HiFiAdapterFilt which is called
# In the data filtering step
git clone https://github.com/sheinasim/HiFiAdapterFilt.git

#Next I ran HiFiAdapter Filt on my raw BAM file from the sequencer with the script below.

#!/bin/bash

#SBATCH --job-name hififilt
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bamtools/2.5.2
module load gcc/11.2.0
module load blast-plus/2.12.0

TOOL_DIR=/anvil/scratch/x-srudnick/Woodrat_WGS/HiFiAdapterFilt

export PATH=$PATH:$TOOL_DIR
export PATH=$PATH:$TOOL_DIR/DB

bash "$TOOL_DIR/hifiadapterfilt.sh" -p m84082_260311_175030_s3.hifi_reads.bc2024 -l 44 -m 97 -t 8

#This ran correctly and provided me with the following files:

m84082_260311_175030_s3.hifi_reads.bc2024.blocklist	
0 Bytes	7/6/2026 9:56:15 PM
file	m84082_260311_175030_s3.hifi_reads.bc2024.contaminant.blastout	
2.22 GB	7/6/2026 9:56:11 PM
file	m84082_260311_175030_s3.hifi_reads.bc2024.filt.fastq.gz	
43.4 GB	7/6/2026 10:05:41 PM
file	m84082_260311_175030_s3.hifi_reads.bc2024.stats	
362 Bytes	7/6/2026 10:07:14 PM

#For downstream assembly I will be using the fastqz.gz (filtered file) our stats results retuneds the following results:

Started on Mon Jul  6 17:43:09 EDT 2026
For the m84082_260311_175030_s3.hifi_reads.bc2024 dataset:
  Removing reads containing adapters a minimum of 44 bp in length and 97% match.

Number of ccs reads: 6052869
Number of adapter contaminated ccs reads: 0 (0% of total)
Number of ccs reads retained: 6052869 (100% of total)

Finished on Mon Jul  6 22:07:14 EDT 2026

#This ensured no remnant adapters are present in the reads. 
#To verify that the total read count was not impacted by the filtering I ran this in the terminal directly
x-srudnick@login03.anvil:[Woodrat_WGS] $ bamtools count -in m84082_260311_175030_s3.hifi_reads.bc2024.bam

#this returned that we had a read count of 6052869 which was the exact number of reads reported in the stats file.

#I was satisfied with HiFiAdapterFiltering at this point.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Next I need to perform coverage validation. I am going to start by converting my fastq file into a fasta file to stay closer to the steps outlined in the CCGP steps
#I used the script below to do this:

#!/bin/bash

#SBATCH --job-name=fastq2fasta
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load seqkit/2.3.1
module load htslib/1.17

seqkit fq2fa m84082_260311_175030_s3.hifi_reads.bc2024.filt.fastq.gz -o m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta

bgzip m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta

#Once this finished I had the zipped, fasta converted file:
m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta.gz

#I wanted to quickly ensure that I did not loose any reads during this conversion so I quickly ran this in the terminal:
  
zgrep -c "^>" m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta.gz

#This returned 6052869 reads in th fasta.gz file, this is good. 

#Next I need to index the fasta with samtools to generate coverage validation. To do so I ran the following script:

#!/bin/bash

#SBATCH --job-name=fastaindex
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load samtools/1.17

samtools faidx m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta.gz

#This indexed our fasta and produced the following files:

m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta.gz.fai	
389 MB	7/8/2026 3:41:09 PM
file	m84082_260311_175030_s3.hifi_reads.bc2024.filt.fasta.gz.gzi	
26.2 MB	7/8/2026 3:41:06 PM

# coverage validation assumes the pacbio adapters have been sequenced. Here I realized that the zero contamination may be due to 
# a database failure. After some investigation I found that the HiFiAdapterFilt does not utilize the most recent adapters and is therefore 
# basically pointless. Although some pipelines do not this filtering I decided it best to temporarily pivot while we await the new db from our
# sequencing provider. 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the mean time I will be running profiling using the VGP pipeline.

https://github.com/VGP/vgp-assembly
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#                                                   DATA PREPERATION

#KMER profiling is the first step. I will begin using Meryl to do so. 
#Meryl uses a fastq file as an input. I want to begin by converting the original BAM to do this rather than the fastq from
# the failed adapterfiltering before

#this is the shell i ran to convert the bam to a fastq

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

# This produced this file

m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz	
35 GB	7/10/2026 4:21:04 PM

#To verify no reads were lost during the conversion I ran: (Directly in the terminal)

module load seqkit
seqkit stats m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#This found:

x-srudnick@login07.anvil:[Meryl] $ seqkit stats m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz
file                                                format  type   num_seqs         sum_len  min_len   avg_len  max_len
m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz  FASTQ   DNA   6,052,869  97,170,520,930      166  16,053.6   60,986


#Good results, we did not lose anything in the conversion. Now I can run Meryl confidently.

#                                                 KMER COUNTING WITH MERYL

#Documentation: https://github.com/marbl/meryl/blob/master/README.md
#Guide: https://meryl.readthedocs.io/en/latest/quick-start.html

#To run Meryl I ran the following shell script below:

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

# Meryl finished running sucessfully and produced a large slurm file which analyzed 31 kmers in my database.
#This was the output file for the run:

woodrat.k31.meryl

#Next I want to validate the database before implementing other Meryl steps. 
# To do so I ran:

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

#This ran correctly and produced the Slurm file (slurm-19150265.out) which contains all the validating statistics of our
#database below:

User guides for each biocontainer module can be found in
https://www.rcac.purdue.edu/knowledge/biocontainers 

Found 1 command tree.
Number of 31-mers that are:
unique             2138735615  (exactly one instance of the kmer is in the input)
# Can indicate: true low-coverage genomic regions, some heterozygous sequence, and error-derived k-mers.
# We will generate a histogram next to validate how these fit into our dataset 
distinct           4645318687  (non-redundant kmer sequences in the input)
# Can indicate: unique genomic sequence,repetitive sequence,heterozygous variants,and some sequencing-error-derived k-mers.
present           96988934860  (...)
# Total present kmers in the dataset
missing   4611686013782069217  (non-redundant kmer sequences not in the input)
#There are 4^31 possible DNA 31-mers (about 4.6×10^18). Your dataset contains only a tiny fraction of all mathematically 
#possible 31-mers, so the overwhelming majority are missing.This statistic is expected and is not used for genome profiling

number of   cumulative   cumulative     presence
distinct     fraction     fraction   in dataset
frequency        kmers     distinct        total       (1e-6)
--------- ------------ ------------ ------------ ------------
  1   2138735615       0.4604       0.0221     0.000010
2    150561103       0.4928       0.0252     0.000021
3     31438974       0.4996       0.0261     0.000031
4     11598095       0.5021       0.0266     0.000041
5      5941453       0.5034       0.0269     0.000052
6      3896521       0.5042       0.0272     0.000062
7      3351355       0.5049       0.0274     0.000072
8      3639696       0.5057       0.0277     0.000082
9      4586870       0.5067       0.0281     0.000093
10      6142827       0.5080       0.0288     0.000103
#I have only attatched the first ten lines here

#To validate these results ad visualize how they interact I then genrated a historgram using the following shell script:

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

#This ran and produced:

woodrat.k31.hist	
477 KB	7/13/2026 11:33:05 AM

#This is the key file genomescope will use to generate estimates of our genome size.
# The histogram saves the kmer results like this:
1	2138735615
2	150561103
3	31438974
4	11598095
5	5941453
6	3896521
7	3351355
8	3639696
9	4586870
10	6142827
# These results show Frequency and Number of distinct Kmers respectivley.

#                                                         GENOMESCOPE

# Documentation https://github.com/tbenavi1/genomescope2.0

#To run Genome scope and analyze the historgram which we just created previously, I elected to use the online interface
# which provides the same results with an easier to use interface.

#First I went to http://genomescope.org/genomescope2.0/

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
#(http://genomescope.org/genomescope2.0/analysis.php?code=Jtbymb38K1zjmDiRXcJv)

# The modeled peak was around 34%

# At this point I have generated all analyses following the VGP pipeline. The next step is to assembly the genome with hifiasm

#                                                 HIFIASM (Contiging solo w HiC)

# AGI has confirmed that the HiFi reads do not contain any adapters and the HiC data is suitible for HiFiasm assembly. Today (7/15) I plan to begin assembly with hifiasm
# My plan is to combine the two HiC runs to produce two combines files (r1 and r2) Then to run hifiasm producing a haplotype resolves assembly. 

#Below is the shell script I submitted to concatenate the two HiC files into one combines R1 and R2 file

#!/bin/bash

#SBATCH --job-name=hicmerge
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --time=8:00:00

cd $SLURM_SUBMIT_DIR

cat WOODRAT604_L_R1.fastq.gz WOODRAT_L_R1.fastq.gz > combined_woodrat_L_R1.fastq.gz

cat WOODRAT604_L_R2.fastq.gz WOODRAT_L_R2.fastq.gz > combined_woodrat_L_R2.fastq.gz

echo "Hi-C merge completed successfully."

# This ran and produces the following files:

combined_woodrat_L_R1.fastq.gz	
55.4 GB	7/15/2026 12:44:52 PM
file	combined_woodrat_L_R2.fastq.gz	
55.6 GB	7/15/2026 12:46:09 PM

# To verify these files contain the same sequence counts as the original files I ran this:

#!/bin/bash

#SBATCH --job-name=seqkit
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128GB
#SBATCH --time=2:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load seqkit/2.3.1

seqkit stats combined_woodrat_L_R1.fastq.gz combined_woodrat_L_R2.fastq.gz

# This produced

User guides for each biocontainer module can be found in
https://www.rcac.purdue.edu/knowledge/biocontainers 
file                            format  type     num_seqs          sum_len  min_len  avg_len  max_len
combined_woodrat_L_R1.fastq.gz  FASTQ   DNA   857,426,762  128,614,014,300      150      150      150
combined_woodrat_L_R2.fastq.gz  FASTQ   DNA   857,426,762  128,614,014,300      150      150      150

#These results were the same number of sequences of both my r1 and r2 files previously which to me indicates we did not lose any data when concatenating the files
# At this point I was ready to run hifiasm

#This was the script I ran for HiFiasm

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

hifiasm -o woodrat_hifiandhic.asm -t32 --h1 combined_woodrat_L_R1.fastq.gz --h2 combined_woodrat_L_R2.fastq.gz m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz

#HiFiasm sucessfully ran in about 6 hours. This is roughly the same amount of time it took to run this with the HiFi data alone. 
# The output files for the run are shown below:

woodrat_hifiandhic.asm.ec.bin	
23 GB	7/15/2026 8:16:34 PM
file	woodrat_hifiandhic.asm.ovlp.source.bin	
13 GB	7/15/2026 8:17:28 PM
file	woodrat_hifiandhic.asm.ovlp.reverse.bin	
9.43 GB	7/15/2026 8:18:06 PM
file	woodrat_hifiandhic.asm.hic.r_utg.gfa	
4.61 GB	7/15/2026 8:24:30 PM
file	woodrat_hifiandhic.asm.hic.r_utg.noseq.gfa	
95.6 MB	7/15/2026 8:24:32 PM
file	woodrat_hifiandhic.asm.hic.r_utg.lowQ.bed	
6.81 MB	7/15/2026 8:25:06 PM
file	woodrat_hifiandhic.asm.hic.p_utg.gfa	
4.57 GB	7/15/2026 8:25:31 PM
file	woodrat_hifiandhic.asm.hic.p_utg.noseq.gfa	
95.2 MB	7/15/2026 8:25:32 PM
file	woodrat_hifiandhic.asm.hic.p_utg.lowQ.bed	
6.73 MB	7/15/2026 8:26:05 PM
file	woodrat_hifiandhic.asm.hic.p_ctg.gfa	
2.76 GB	7/15/2026 8:30:01 PM
file	woodrat_hifiandhic.asm.hic.p_ctg.noseq.gfa	
63.6 MB	7/15/2026 8:30:02 PM
file	woodrat_hifiandhic.asm.hic.p_ctg.lowQ.bed	
1.56 MB	7/15/2026 8:30:22 PM
file	woodrat_hifiandhic.asm.hic.tlb.bin	
74.6 GB	7/15/2026 8:39:11 PM
file	woodrat_hifiandhic.asm.hic.lk.bin	
5.44 GB	7/15/2026 9:18:30 PM
file	woodrat_hifiandhic.asm.hic.hap1.p_ctg.gfa	
2.5 GB	7/15/2026 9:22:01 PM
file	woodrat_hifiandhic.asm.hic.hap1.p_ctg.noseq.gfa	
57.1 MB	7/15/2026 9:22:02 PM
file	woodrat_hifiandhic.asm.hic.hap1.p_ctg.lowQ.bed	
1.88 MB	7/15/2026 9:22:21 PM
file	woodrat_hifiandhic.asm.hic.hap2.p_ctg.gfa	
2.79 GB	7/15/2026 9:22:31 PM
file	woodrat_hifiandhic.asm.hic.hap2.p_ctg.noseq.gfa	
63.4 MB	7/15/2026 9:22:32 PM
file	woodrat_hifiandhic.asm.hic.hap2.p_ctg.lowQ.bed	
1.78 MB	7/15/2026 9:22:53 PM
file	slurm-19255719.out	
36.8 KB	7/15/2026 9:23:09 PM

# This all looks good for output and I do not observe any error messages. Now we need to run QC to ensure there were no
# hidden errors. To Qc our hap1 and hap 2 assembly we will use gfastats. 

#                                                   GFAstats

#gfastats is a tool suite that allows for manipulation and evaluation of FASTA and GFA files, but in this instance we 
# will use it to convert our GFAs to FASTA files. Later on, we will use it to generate standard summary statistics for 
#our assemblies.

#Documentation: https://github.com/vgl-hub/gfastats

# To run gfastats I first created a new directory with the two phased haplotype assemblies as was done in the VGP workflow 
# to run WF4.

#These are the two files
woodrat_hifiandhic.asm.hic.hap1.p_ctg.gfa	
2.5 GB	7/16/2026 12:45:43 PM
file	woodrat_hifiandhic.asm.hic.hap2.p_ctg.gfa	
2.79 GB	7/16/2026 12:46:55 PM

#We will run gfa in two steps as is done in the VGP pipeline. First is the fasta conversion. 
#Below is the shell script i ran for this first conversion

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

# This produced the fasta converted files for hap 1 and hap 2 below:

woodrat_hap1.fasta	
2.44 GB	7/16/2026 1:29:31 PM

woodrat_hap2.fasta	
2.73 GB	7/16/2026 1:29:49 PM #Haplotype 2 is larger?

slurm-19298174.out	
2.31 KB	7/16/2026 1:29:49 PM

# These files were interesting in size but i continued to stats. I then ran the following scripts to generate statistics 

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

# This ran sucessfully and produced the following files:
slurm-19298591.out	
110 Bytes	7/16/2026 1:44:50 PM

file	woodrat_hap1_gfastats.txt	
863 Bytes	7/16/2026 1:45:06 PM

file	woodrat_hap2_gfastats.txt	
863 Bytes	7/16/2026 1:45:22 PM

#Haplotype 1 summary results:

+++Assembly summary+++: 
  # scaffolds: 312
  Total scaffold length: 2624980539
Average scaffold length: 8413399.16
Scaffold N50: 35249638
Scaffold auN: 39223716.75
Scaffold L50: 25
Largest scaffold: 94591519
Smallest scaffold: 20092
# contigs: 312
Total contig length: 2624980539
Average contig length: 8413399.16
Contig N50: 35249638
Contig auN: 39223716.75
Contig L50: 25
Largest contig: 94591519
Smallest contig: 20092
# gaps in scaffolds: 0
Total gap length in scaffolds: 0
Average gap length in scaffolds: 0.00
Gap N50 in scaffolds: 0
Gap auN in scaffolds: 0.00
Gap L50 in scaffolds: 0
Largest gap in scaffolds: 0
Smallest gap in scaffolds: 0
Base composition (A:C:G:T): 771077500:541084690:541538105:771280244
GC content %: 41.24
# soft-masked bases: 0
# segments: 312
Total segment length: 2624980539
Average segment length: 8413399.16
# gaps: 0
# paths: 312

#Haplotype 2 summary:

+++Assembly summary+++: 
  # scaffolds: 299
  Total scaffold length: 2931396252
Average scaffold length: 9804000.84
Scaffold N50: 35204138
Scaffold auN: 39970146.86
Scaffold L50: 27
Largest scaffold: 95463868
Smallest scaffold: 26626
# contigs: 299
Total contig length: 2931396252
Average contig length: 9804000.84
Contig N50: 35204138
Contig auN: 39970146.86
Contig L50: 27
Largest contig: 95463868
Smallest contig: 26626
# gaps in scaffolds: 0
Total gap length in scaffolds: 0
Average gap length in scaffolds: 0.00
Gap N50 in scaffolds: 0
Gap auN in scaffolds: 0.00
Gap L50 in scaffolds: 0
Largest gap in scaffolds: 0
Smallest gap in scaffolds: 0
Base composition (A:C:G:T): 864935993:600884121:600985274:864590864
GC content %: 41.00
# soft-masked bases: 0
# segments: 299
Total segment length: 2931396252
Average segment length: 9804000.84
# gaps: 0
# paths: 299

# Next we rerun gfa stats with the previously generated estimated genome size from genome scope. This will allow us to generate
# ng50 stats. To do this we rerun the gfastats stat generatio with an additional line 
# Genomescope outputs a range of a min and a max genome size. The yreported using the max genome size estimate for their input
# My maximum genome size was 2,743,204,927 bp 

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
2743204927 \
--stats \
> woodrat_hap1_gfastats_egs.txt

gfastats \
-f woodrat_hap2.fasta \
2743204927 \
--stats \
> woodrat_hap2_gfastats_egs.txt

# After running the stats again with estimated genome size, these were the expanded results:

#Results Hap 1:

+++Assembly summary+++: 
  Expected genome size: 2743204927
# scaffolds: 312
Total scaffold length: 2624980539
Average scaffold length: 8413399.16
Scaffold N50: 35249638
Scaffold auN: 39223716.75
Scaffold L50: 25
Scaffold NG50: 32780092
Scaffold auNG: 37533285.29
Scaffold LG50: 27
Largest scaffold: 94591519
Smallest scaffold: 20092
# contigs: 312
Total contig length: 2624980539
Average contig length: 8413399.16
Contig N50: 35249638
Contig auN: 39223716.75
Contig L50: 25
Contig NG50: 32780092
Contig auNG: 37533285.29
Contig LG50: 27
Largest contig: 94591519
Smallest contig: 20092
# gaps in scaffolds: 0
Total gap length in scaffolds: 0
Average gap length in scaffolds: 0.00
Gap N50 in scaffolds: 0
Gap auN in scaffolds: 0.00
Gap L50 in scaffolds: 0
Largest gap in scaffolds: 0
Smallest gap in scaffolds: 0
Base composition (A:C:G:T): 771077500:541084690:541538105:771280244
GC content %: 41.24
# soft-masked bases: 0
# segments: 312
Total segment length: 2624980539
Average segment length: 8413399.16
# gaps: 0
# paths: 312


#Results Hap 2:

+++Assembly summary+++: 
  Expected genome size: 2743204927
# scaffolds: 299
Total scaffold length: 2931396252
Average scaffold length: 9804000.84
Scaffold N50: 35204138
Scaffold auN: 39970146.86
Scaffold L50: 27
Scaffold NG50: 35460640
Scaffold auNG: 42712207.73
Scaffold LG50: 25
Largest scaffold: 95463868
Smallest scaffold: 26626
# contigs: 299
Total contig length: 2931396252
Average contig length: 9804000.84
Contig N50: 35204138
Contig auN: 39970146.86
Contig L50: 27
Contig NG50: 35460640
Contig auNG: 42712207.73
Contig LG50: 25
Largest contig: 95463868
Smallest contig: 26626
# gaps in scaffolds: 0
Total gap length in scaffolds: 0
Average gap length in scaffolds: 0.00
Gap N50 in scaffolds: 0
Gap auN in scaffolds: 0.00
Gap L50 in scaffolds: 0
Largest gap in scaffolds: 0
Smallest gap in scaffolds: 0
Base composition (A:C:G:T): 864935993:600884121:600985274:864590864
GC content %: 41.00
# soft-masked bases: 0
# segments: 299
Total segment length: 2931396252
Average segment length: 9804000.84
# gaps: 0
# paths: 299

# I followed the VGP pipeline by joining the two summary files. To do so I ran the following in terminal:

paste woodrat_hap1_gfastats_egs.txt woodrat_hap2_gfastats_egs.txt > woodrat_gfastats_hap1_hap2_full.txt

# I now have a joined summary file. 
# The VGP pipeline asks: 
#What is the length of the longest contigs in the assemblies? 
Hap1 Largest contig: 94591519	Hap2 Largest contig: 95463868

#What are the N50 values of the two assemblies? Are they very different from each other?
Hap1 Contig N50: 35249638	Hap2 Contig N50: 35204138


#                                                       BUSCO
#Documentation: https://github.com/metashot/busco

#Now I ran busco to identiyf the benchmarking single copy orthologs present in my assembly 

#This was the script 

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
module load busco/5.8.2-py

busco \
-i woodrat_hap1.fasta \
-m genome \
-l mammalia_odb10 \
--metaeuk \
-o busco_hap1 \
-c 32

busco \
-i woodrat_hap2.fasta \
-m genome \
-l mammalia_odb10 \
--metaeuk \
-o busco_hap2 \
-c 32

# The BUSCO analysis revealed that Hap 2 conserved more genes than Hap2. Both assemblies show high completedness greater than 94%, low fragmentation (<1%) and low duplication
# <2%, consistent with my other QC steps. Below are the output files:

Hap 1: 
  logs	
-	7/17/2026 5:39:00 PM
dir	run_mammalia_odb10	
-	7/17/2026 5:39:00 PM
file	short_summary.specific.mammalia_odb10.busco_hap1.txt	
963 Bytes	7/17/2026 5:39:00 PM
file	short_summary.specific.mammalia_odb10.busco_hap1.json

Hap 2:
  logs	
-	7/17/2026 10:08:30 PM
dir	run_mammalia_odb10	
-	7/17/2026 10:08:30 PM
file	short_summary.specific.mammalia_odb10.busco_hap2.json	
2.88 KB	7/17/2026 10:08:30 PM
file	short_summary.specific.mammalia_odb10.busco_hap2.txt	
963 Bytes	7/17/2026 10:08:30 PM

# After running BUSCO the next step is to run QC with Merqury

#                                                           MERQURY
#Documentation: https://github.com/marbl/merqury


#Below is the shell script I used to run merqury. Merqury uses the kmer db generated earlier during meryl


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

# merqury ran in about 20 minutes and produced the following files: 
logs-	7/17/2026 1:29:15 PM
file	slurm-19322333.out	
307 Bytes	7/17/2026 1:29:15 PM
file	woodrat.k31.hist	
477 KB	7/17/2026 1:29:15 PM
file	woodrat.k31.hist.ploidy	
44 Bytes	7/17/2026 1:29:16 PM
file	woodrat.k31.filt	
2 Bytes	7/17/2026 1:29:38 PM
dir	woodrat_hap1.meryl	
-	7/17/2026 1:30:22 PM
file	woodrat_merqury.woodrat_hap1.spectra-cn.hist	
725 KB	7/17/2026 1:32:47 PM
file	woodrat_merqury.woodrat_hap1.only.hist	
17 Bytes	7/17/2026 1:33:08 PM
file	woodrat_merqury.woodrat_hap1.spectra-cn.ln.png	
108 KB	7/17/2026 1:33:10 PM
file	woodrat_merqury.woodrat_hap1.spectra-cn.fl.png	
107 KB	7/17/2026 1:33:11 PM
file	woodrat_merqury.woodrat_hap1.spectra-cn.st.png	
121 KB	7/17/2026 1:33:34 PM
file	woodrat_merqury.woodrat_hap1.qv	
11.7 KB	7/17/2026 1:34:05 PM
file	woodrat_hap1_only.bed	
450 KB	7/17/2026 1:34:53 PM
file	woodrat_hap1_only.wig	
380 KB	7/17/2026 1:35:22 PM
dir	woodrat_hap2.meryl	
-	7/17/2026 1:36:10 PM
file	woodrat_merqury.woodrat_hap2.spectra-cn.hist	
715 KB	7/17/2026 1:38:54 PM
file	woodrat_merqury.woodrat_hap2.only.hist	
16 Bytes	7/17/2026 1:39:18 PM
file	woodrat_merqury.woodrat_hap2.spectra-cn.fl.png	
104 KB	7/17/2026 1:39:20 PM
file	woodrat_merqury.woodrat_hap2.spectra-cn.ln.png	
104 KB	7/17/2026 1:39:20 PM
file	woodrat_merqury.woodrat_hap2.spectra-cn.st.png	
120 KB	7/17/2026 1:39:43 PM
file	woodrat_merqury.woodrat_hap2.qv	
11.4 KB	7/17/2026 1:40:08 PM
file	woodrat_hap2_only.bed	
408 KB	7/17/2026 1:40:52 PM
file	woodrat_hap2_only.wig	
363 KB	7/17/2026 1:41:18 PM
file	woodrat_merqury.completeness.stats	
133 Bytes	7/17/2026 1:41:57 PM
file	woodrat_merqury.spectra-cn.hist	
708 KB	7/17/2026 1:44:27 PM
file	woodrat_merqury.only.hist	
19 Bytes	7/17/2026 1:44:49 PM
file	woodrat_merqury.dist_only.hist	
64 Bytes	7/17/2026 1:44:49 PM
file	woodrat_merqury.spectra-cn.ln.png	
110 KB	7/17/2026 1:44:51 PM
file	woodrat_merqury.spectra-cn.fl.png	
109 KB	7/17/2026 1:44:51 PM
file	woodrat_merqury.qv	
141 Bytes	7/17/2026 1:45:14 PM
file	woodrat_merqury.spectra-cn.st.png	
124 KB	7/17/2026 1:45:14 PM
file	woodrat_merqury.spectra-asm.hist	
955 KB	7/17/2026 1:46:48 PM
file	woodrat_merqury.spectra-asm.fl.png	
122 KB	7/17/2026 1:46:50 PM
file	woodrat_merqury.spectra-asm.ln.png	
124 KB	7/17/2026 1:46:50 PM
file	woodrat_merqury.spectra-asm.st.png	
129 KB	7/17/2026 1:47:13 PM

# Here are some of the most important reports from Merqury:

#Merqury stats
#From: woodrat_merqury.completeness.stats
#name #all trusted kmers #kmers recovered #expected kmer count #completedness
woodrat_hap1	all	2129667657	2299795571	92.6025
woodrat_hap2	all	2253189993	2299795571	97.9735 # <----------Completeness score
both	all	2292831434	2299795571	99.6972

#Merqury plots
#(.st	standard).(fl	full length,Shows the full k-mer spectrum including the long tail) (.ln	log-normalized	Uses log scaling)
Main plot: 	woodrat_merqury.spectra-asm.st.png #Assembly spectra 
"Merqury ASM plot. This plot tracks the multiplicity of each k-mer found in the HiFi read set and colors it according to which
assemblies contain those k-mers. This can tell you which k-mers are found in only one assembly or shared between them."

#Merqury consensus accuracy quality value (QV) stats. 
# From: 	woodrat_merqury.qv
#Name #error kmers #assembly size #QV (Phred consensus accuracy) #estimated base error
woodrat_hap1	13640	2624971179	67.7567	1.67621e-07
woodrat_hap2	12465	2931387282	68.6274	1.3717e-07 #Very very strong results
Both	26105	5556358461	68.1943	1.51556e-07

#                                               PURGING

#Documentation: https://github.com/dfguan/purge_dups

#The Vertebrate genome project typically does not necessitate purging but my genome is a good cantidate to do this just to be extra sure based on the duplicated material 
# Identified in both BUSCO and merqury. We will be doing this to test whether the phased assemblies contain redundant 
# haplotypic sequence

# The first step was to create a new wd called purge_dups. I first brought in these files 

m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz	<-- This was made 
35 GB	7/27/2026 3:24:44 PM
file	woodrat_hap1.fasta	
2.44 GB	7/27/2026 3:28:29 PM
file	woodrat_hap2.fasta	
2.73 GB	7/27/2026 3:28:36 PM

# Step 1A: Read-depth calculation: 
#The first step is to allign the hap1 fasta to the raw hifi reads and then repeat this for hap2:

# here is the script I ran to do so:

#!/bin/bash

#SBATCH --job-name=allign_h1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28

minimap2 -xmap-hifi -I 8G woodrat_hap1.fasta m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz | gzip -c - > hap1.paf.gz #IF I RERUN REMEMBER THAT I CHANGED -X MAP TO -XMAP LIKE IN GITHUB

# I also ran this for hap 2 at the same time:

#!/bin/bash

#SBATCH --job-name=allign_h2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28

minimap2 -x map-hifi -I 8G woodrat_hap2.fasta m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz | gzip -c - > hap2.paf.gz

# This finished in roughly 4 hours. Hap2 took much longer and finished in 4.5 hours

#This was the hap 1 output:

hap1.paf.gz	
387 MB	7/28/2026 2:15:37 PM
file	slurm-19550659.out	
12 KB	7/28/2026 2:15:37 PM


#This was the hap 2 output:

hap2.paf.gz	
345 MB	7/28/2026 3:34:11 PM
file	slurm-19550665.out	
12.1 KB	7/28/2026 3:34:11 PM

#The next step to run purge_dups is to produce pbcstat to do so I ran the pbc steps with the shells below for h1 and h2:

#!/bin/bash

#SBATCH --job-name=pbc_h1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load purge_dups/1.2.6_fixed

pbcstat hap1.paf.gz

#!/bin/bash

#SBATCH --job-name=pbc_h2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load purge_dups/1.2.6_fixed

pbcstat hap2.paf.gz


# I ended up having to create a dedicated subdirectory for each haplotype since this output for pbc has a fixed output file 
# that will overwrite the other 

#This did sucessfully produce the same three files for each haplotype respectivley:

#Hap1: 
PB.base.cov	
223 MB	7/28/2026 4:39:52 PM
file	PB.cov.wig	
7.37 MB	7/28/2026 4:39:52 PM
file	PB.stat	
4.5 KB	7/28/2026 4:39:50 PM
file	slurm-19555647.out	
155 Bytes	7/28/2026 4:39:52 PM

#Hap2: 
PB.base.cov	
225 MB	7/28/2026 4:39:49 PM
file	PB.cov.wig	
8.2 MB	7/28/2026 4:39:49 PM
file	PB.stat	
3.85 KB	7/28/2026 4:39:48 PM
file	slurm-19555646.out	
155 Bytes	7/28/2026 4:39:49 PM

# I then ran the calcuts function for each haplotype with this script:

#!/bin/bash

#SBATCH --job-name=calcut_hap2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load purge_dups/1.2.6_fixed

calcuts PB.stat > cutoffs 2>calcults.log

#!/bin/bash

#SBATCH --job-name=calcut_hap1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load purge_dups/1.2.6_fixed

calcuts PB.stat > cutoffs 2>calcults.log

#This produced:

calcults.log	
199 Bytes	7/28/2026 4:57:07 PM
file	cutoffs	
17 Bytes	7/28/2026 4:57:07 PM
file	slurm-19555696.out	
110 Bytes	7/28/2026 4:57:07 PM

#and for hap1: 

slurm-19555695.out	
110 Bytes	7/28/2026 4:57:08 PM
file	cutoffs	
18 Bytes	7/28/2026 4:57:09 PM
file	calcults.log	
195 Bytes	7/28/2026 4:57:09 PM

#Step 1B. Now that we have completed 1A we move to Spliting the assembly and self-alignment:

#!/bin/bash

#SBATCH --job-name split_h1
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=16:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

split_fa woodrat_hap1.fasta > woodrat_hap1.fasta.split

minimap2 -xasm5 -DP woodrat_hap1.fasta.split woodrat_hap1.fasta.split | gzip -c > hap1.self.paf.gz

#!/bin/bash

#SBATCH --job-name split_h2
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=16:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

split_fa woodrat_hap2.fasta > woodrat_hap2.fasta.split

minimap2 -xasm5 -DP woodrat_hap2.fasta.split woodrat_hap2.fasta.split | gzip -c > hap2.self.paf.gz

# This produced: 
hap1.self.paf.gz	
856 MB	7/29/2026 11:40:10 AM

hap2.self.paf.gz	
1.38 GB	7/29/2026 11:57:08 AM

#Step 2 is to now purge the duplicates. To do so I ran these scripts:

#!/bin/bash

#SBATCH --job-name purge_h1
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=16:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

purge_dups -2 -T h1_cutoffs -c PB.h1base.cov hap1.self.paf.gz > dups_h1.bed 2> purge_h1_dups.log

#!/bin/bash

#SBATCH --job-name purge_h2
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=16:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

purge_dups -2 -T h2_cutoffs -c PB.h2base.cov hap2.self.paf.gz > dups_h2.bed 2> purge_h2_dups.log

# This finished very fast and produced 
slurm-19565186.out	
110 Bytes	7/29/2026 1:38:32 PM
file	dups_h1.bed	
12.7 KB	7/29/2026 1:39:03 PM
file	purge_h1_dups.log	
653 KB	7/29/2026 1:39:03 PM

#and for h2:
slurm-19565183.out	
110 Bytes	7/29/2026 1:38:08 PM
file	purge_h2_dups.log	
1.42 MB	7/29/2026 1:38:57 PM
file	dups_h2.bed	
9.08 KB	7/29/2026 1:38:58 PM

# Step 3: Get purged primary and haplotig sequences from draft assembly.

#!/bin/bash

#SBATCH --job-name getseq_h1
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

get_seqs -e dups_h1.bed woodrat_hap1.fasta


#!/bin/bash

#SBATCH --job-name getseq_h2
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB
#SBATCH --time=1:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load minimap2/2.28
module load purge_dups/1.2.6_fixed

get_seqs -e dups_h2.bed woodrat_hap2.fasta

#This also finished very quickly. Now we have two purged assemblies: 

file	woodrat_hap1_purged.fa	
2.25 GB	7/29/2026 1:45:38 PM
file	hap.fa	
197 MB	7/29/2026 1:45:38 PM
file	slurm-19565221.out	
3.56 KB	7/29/2026 1:45:38 PM

slurm-19565224.out	
864 Bytes	7/29/2026 1:46:37 PM
file	woodrat_hap2_purged.fa	
2.47 GB	7/29/2026 1:46:37 PM
file	hap.fa	
268 MB	7/29/2026 1:46:37 PM

# Purge_dups is finished now. Now we must validate the purged fasta file with gfa stats, merqury and busco. Since we have
# basically already done this now we must compare the results and see which performs better to use for scaffolding:

# I'm going to do these steps 1 haplotype at a time. Starting with Hap 1

#                                                                   GFASTATS


#!/bin/bash

#SBATCH --job-name=purgegfa
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
-f woodrat_hap1_purged.fa \
--stats \
> purged_hap1_gfastats.txt

gfastats \
-f woodrat_hap2_purged.fa \
--stats \
> purged_hap2_gfastats.txt

#This produced:
slurm-19579863.out	
110 Bytes	7/30/2026 1:49:05 PM
file	purged_hap1_gfastats.txt	
866 Bytes	7/30/2026 1:49:20 PM
file	purged_hap2_gfastats.txt	
866 Bytes	7/30/2026 1:49:35 PM

# The output files reported okay results. We lost a lot of assembvly size but had fewer contigs and beter N50 value. Now we will perform BUSCO
# To identify if the loss of size reduces out overall completeness valeus 

#                                                                         Purged_BUSCO 

#This is the shell script i ran:

#!/bin/bash

#SBATCH --job-name=purge_BUSCO
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load busco/5.8.2-py

busco \
-i woodrat_hap1_purged.fa \
-m genome \
-l mammalia_odb10 \
--metaeuk \
-o busco_purg_hap1 \
-c 32

busco \
-i woodrat_hap2_purged.fa \
-m genome \
-l mammalia_odb10 \
--metaeuk \
-o busco_purg_hap2 \
-c 32

# The output SLURM file from this run summarizes how the BUSCO went when analyzing the purged results:

User guides for each biocontainer module can be found in
https://www.rcac.purdue.edu/knowledge/biocontainers 
Please remember to set AUGUSTUS_CONFIG_PATH - see note in 'module help
busco/5.8.2-py' 
2026-07-30 13:59:13 INFO:	***** Start a BUSCO v5.8.2 analysis, current time: 07/30/2026 13:59:13 *****
  2026-07-30 13:59:13 INFO:	Configuring BUSCO with local environment
2026-07-30 13:59:13 INFO:	Running genome mode
2026-07-30 13:59:13 INFO:	Downloading information on latest versions of BUSCO data...
2026-07-30 13:59:14 INFO:	Input file is /anvil/scratch/x-srudnick/Woodrat_WGS/purge_dups/purged_BUSCO/woodrat_hap1_purged.fa
2026-07-30 13:59:14 INFO:	Downloading file 'https://busco-data.ezlab.org/v5/data/lineages/mammalia_odb10.2024-01-08.tar.gz'
2026-07-30 13:59:27 INFO:	Decompressing file '/anvil/scratch/x-srudnick/Woodrat_WGS/purge_dups/purged_BUSCO/busco_downloads/lineages/mammalia_odb10.tar.gz'
2026-07-30 14:00:08 INFO:	Running BUSCO using lineage dataset mammalia_odb10 (eukaryota, 2024-01-08)
2026-07-30 14:00:08 INFO:	Running 1 job(s) on bbtools, starting at 07/30/2026 14:00:08
2026-07-30 14:00:23 INFO:	[bbtools]	1 of 1 task(s) completed
2026-07-30 14:00:24 INFO:	Running 1 job(s) on metaeuk, starting at 07/30/2026 14:00:24
2026-07-30 15:56:16 INFO:	[metaeuk]	1 of 1 task(s) completed
2026-07-30 15:56:18 INFO:	***** Run HMMER on gene sequences *****
  2026-07-30 15:56:18 INFO:	Running 9226 job(s) on hmmsearch, starting at 07/30/2026 15:56:18
2026-07-30 15:56:36 INFO:	[hmmsearch]	923 of 9226 task(s) completed
2026-07-30 15:56:51 INFO:	[hmmsearch]	1846 of 9226 task(s) completed
2026-07-30 15:57:04 INFO:	[hmmsearch]	2768 of 9226 task(s) completed
2026-07-30 15:57:14 INFO:	[hmmsearch]	3691 of 9226 task(s) completed
2026-07-30 15:57:23 INFO:	[hmmsearch]	4614 of 9226 task(s) completed
2026-07-30 15:57:36 INFO:	[hmmsearch]	5536 of 9226 task(s) completed
2026-07-30 15:58:00 INFO:	[hmmsearch]	6459 of 9226 task(s) completed
2026-07-30 15:58:26 INFO:	[hmmsearch]	7381 of 9226 task(s) completed
2026-07-30 15:58:46 INFO:	[hmmsearch]	8304 of 9226 task(s) completed
2026-07-30 15:59:23 INFO:	[hmmsearch]	9226 of 9226 task(s) completed
2026-07-30 16:00:06 INFO:	734 candidate overlapping regions found
2026-07-30 16:00:06 INFO:	108541 exons in total
2026-07-30 16:02:08 INFO:	Extracting missing and fragmented buscos from the file refseq_db.faa...
2026-07-30 16:02:40 INFO:	Running 1 job(s) on metaeuk, starting at 07/30/2026 16:02:40
2026-07-30 18:20:33 INFO:	[metaeuk]	1 of 1 task(s) completed
2026-07-30 18:20:34 INFO:	***** Run HMMER on gene sequences *****
  2026-07-30 18:20:34 INFO:	Running 665 job(s) on hmmsearch, starting at 07/30/2026 18:20:34
2026-07-30 18:20:43 INFO:	[hmmsearch]	67 of 665 task(s) completed
2026-07-30 18:20:49 INFO:	[hmmsearch]	133 of 665 task(s) completed
2026-07-30 18:20:56 INFO:	[hmmsearch]	200 of 665 task(s) completed
2026-07-30 18:21:03 INFO:	[hmmsearch]	266 of 665 task(s) completed
2026-07-30 18:21:09 INFO:	[hmmsearch]	333 of 665 task(s) completed
2026-07-30 18:21:16 INFO:	[hmmsearch]	399 of 665 task(s) completed
2026-07-30 18:21:23 INFO:	[hmmsearch]	466 of 665 task(s) completed
2026-07-30 18:21:30 INFO:	[hmmsearch]	532 of 665 task(s) completed
2026-07-30 18:21:37 INFO:	[hmmsearch]	599 of 665 task(s) completed
2026-07-30 18:21:43 INFO:	[hmmsearch]	665 of 665 task(s) completed
2026-07-30 18:21:49 INFO:	743 candidate overlapping regions found
2026-07-30 18:21:49 INFO:	95002 exons in total
2026-07-30 18:23:38 INFO:	Results:	C:94.3%[S:92.9%,D:1.4%],F:0.7%,M:5.0%,n:9226	   

2026-07-30 18:24:12 INFO:	
  
  ---------------------------------------------------
  |Results from dataset mammalia_odb10               |
  ---------------------------------------------------
  |C:94.3%[S:92.9%,D:1.4%],F:0.7%,M:5.0%,n:9226      |
  |8704    Complete BUSCOs (C)                       |
  |8572    Complete and single-copy BUSCOs (S)       |
  |132    Complete and duplicated BUSCOs (D)         |
  |63    Fragmented BUSCOs (F)                       |
  |459    Missing BUSCOs (M)                         |
  |9226    Total BUSCO groups searched               |
  ---------------------------------------------------
  2026-07-30 18:24:12 INFO:	BUSCO analysis done. Total running time: 15897 seconds
2026-07-30 18:24:12 INFO:	Results written in /anvil/scratch/x-srudnick/Woodrat_WGS/purge_dups/purged_BUSCO/busco_purg_hap1
2026-07-30 18:24:12 INFO:	For assistance with interpreting the results, please consult the userguide: https://busco.ezlab.org/busco_userguide.html

2026-07-30 18:24:12 INFO:	Visit this page https://gitlab.com/ezlab/busco#how-to-cite-busco to see how to cite BUSCO
2026-07-30 18:24:18 INFO:	***** Start a BUSCO v5.8.2 analysis, current time: 07/30/2026 18:24:18 *****
  2026-07-30 18:24:18 INFO:	Configuring BUSCO with local environment
2026-07-30 18:24:18 INFO:	Running genome mode
2026-07-30 18:24:18 INFO:	Downloading information on latest versions of BUSCO data...
2026-07-30 18:24:20 INFO:	Input file is /anvil/scratch/x-srudnick/Woodrat_WGS/purge_dups/purged_BUSCO/woodrat_hap2_purged.fa
2026-07-30 18:24:23 INFO:	Running BUSCO using lineage dataset mammalia_odb10 (eukaryota, 2024-01-08)
2026-07-30 18:24:23 INFO:	Running 1 job(s) on bbtools, starting at 07/30/2026 18:24:23
2026-07-30 18:24:38 INFO:	[bbtools]	1 of 1 task(s) completed
2026-07-30 18:24:38 INFO:	Running 1 job(s) on metaeuk, starting at 07/30/2026 18:24:38
2026-07-30 20:24:19 INFO:	[metaeuk]	1 of 1 task(s) completed
2026-07-30 20:24:21 INFO:	***** Run HMMER on gene sequences *****
  2026-07-30 20:24:21 INFO:	Running 9226 job(s) on hmmsearch, starting at 07/30/2026 20:24:21
2026-07-30 20:24:41 INFO:	[hmmsearch]	923 of 9226 task(s) completed
2026-07-30 20:24:58 INFO:	[hmmsearch]	1846 of 9226 task(s) completed
2026-07-30 20:25:12 INFO:	[hmmsearch]	2768 of 9226 task(s) completed
2026-07-30 20:25:24 INFO:	[hmmsearch]	3691 of 9226 task(s) completed
2026-07-30 20:25:34 INFO:	[hmmsearch]	4614 of 9226 task(s) completed
2026-07-30 20:25:49 INFO:	[hmmsearch]	5536 of 9226 task(s) completed
2026-07-30 20:26:15 INFO:	[hmmsearch]	6459 of 9226 task(s) completed
2026-07-30 20:26:44 INFO:	[hmmsearch]	7381 of 9226 task(s) completed
2026-07-30 20:27:07 INFO:	[hmmsearch]	8304 of 9226 task(s) completed
2026-07-30 20:27:46 INFO:	[hmmsearch]	9226 of 9226 task(s) completed
2026-07-30 20:28:32 INFO:	619 candidate overlapping regions found
2026-07-30 20:28:32 INFO:	111432 exons in total
2026-07-30 20:30:28 INFO:	Extracting missing and fragmented buscos from the file refseq_db.faa...
2026-07-30 20:30:51 INFO:	Running 1 job(s) on metaeuk, starting at 07/30/2026 20:30:51
2026-07-30 22:27:59 INFO:	[metaeuk]	1 of 1 task(s) completed
2026-07-30 22:27:59 INFO:	***** Run HMMER on gene sequences *****
  2026-07-30 22:27:59 INFO:	Running 455 job(s) on hmmsearch, starting at 07/30/2026 22:27:59
2026-07-30 22:28:09 INFO:	[hmmsearch]	46 of 455 task(s) completed
2026-07-30 22:28:16 INFO:	[hmmsearch]	91 of 455 task(s) completed
2026-07-30 22:28:24 INFO:	[hmmsearch]	137 of 455 task(s) completed
2026-07-30 22:28:31 INFO:	[hmmsearch]	182 of 455 task(s) completed
2026-07-30 22:28:39 INFO:	[hmmsearch]	228 of 455 task(s) completed
2026-07-30 22:28:47 INFO:	[hmmsearch]	273 of 455 task(s) completed
2026-07-30 22:28:54 INFO:	[hmmsearch]	319 of 455 task(s) completed
2026-07-30 22:29:01 INFO:	[hmmsearch]	364 of 455 task(s) completed
2026-07-30 22:29:09 INFO:	[hmmsearch]	410 of 455 task(s) completed
2026-07-30 22:29:18 INFO:	[hmmsearch]	455 of 455 task(s) completed
2026-07-30 22:29:22 INFO:	809 candidate overlapping regions found
2026-07-30 22:29:22 INFO:	97110 exons in total
2026-07-30 22:31:24 INFO:	Results:	C:96.7%[S:95.1%,D:1.6%],F:0.7%,M:2.7%,n:9226	   

2026-07-30 22:31:54 INFO:	
  
  ---------------------------------------------------
  |Results from dataset mammalia_odb10               |
  ---------------------------------------------------
  |C:96.7%[S:95.1%,D:1.6%],F:0.7%,M:2.7%,n:9226      |
  |8919    Complete BUSCOs (C)                       |
  |8775    Complete and single-copy BUSCOs (S)       |
  |144    Complete and duplicated BUSCOs (D)         |
  |61    Fragmented BUSCOs (F)                       |
  |246    Missing BUSCOs (M)                         |
  |9226    Total BUSCO groups searched               |
  ---------------------------------------------------
  2026-07-30 22:31:54 INFO:	BUSCO analysis done. Total running time: 14854 seconds
2026-07-30 22:31:54 INFO:	Results written in /anvil/scratch/x-srudnick/Woodrat_WGS/purge_dups/purged_BUSCO/busco_purg_hap2
2026-07-30 22:31:54 INFO:	For assistance with interpreting the results, please consult the userguide: https://busco.ezlab.org/busco_userguide.html

2026-07-30 22:31:54 INFO:	Visit this page https://gitlab.com/ezlab/busco#how-to-cite-busco to see how to cite BUSCO

# We see an improvment in almost every busco measurment. The most notable change from our original files to our purged assemblies
# is the number of contigs. Hap1: 312 → 122 contigs Hap2: 350 → 144 contigs. At this stage I need to decide if these purged files should be used
# instread of the original files for scaffolding.


# The next validation will be run with merqury

#                                                       MERQURY_PURGED

#Below is the shell script I used to run merqury. Merqury uses the kmer db generated earlier during meryl

#Here is the same script I ran peviously adapted to the purged assembly files:


#!/bin/bash

#SBATCH --job-name=purged_Merqury
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
woodrat_hap1_purged.fa \
woodrat_hap2_purged.fa \
woodrat_purged_merqury

#This ran in roughly 15 minutes and produces all cooresponding files as the first merqury run 

# Ultimatley we saw 
# Before purge: Larger assemblies, More contigs, Slightly higher k-mer completeness #More duplicated sequence likely
# After purge:
 312 → 122 contigs (Hap1)
~350 → 144 contigs (Hap2)
N50 slightly improved
BUSCO maintained/improved:
Hap1: 94.3%
Hap2: 96.7%
Merqury QV essentially unchanged
Merqury completeness still extremely high: Combined: 98.85%

# This is good justification to use the purged fastas for scaffolding


# Following th vgp pipeline we should now be at the point to begin scaffolding. With our data types we should be performing HiC scaffolding with YAHS

#                                                                             SCAFFOLDING


#---------------------------------------------------------------------------------------------------------------------------------------------

# Pre-processing Hi-C data. For Scaffolding I will use the purged data. I will run each haplotype independently. 


#                                                                             BWA MEM

#Step 1: Mapping Hi-C reads. To do this I used BWA MEM. bwa mem2 is not available in anvil

#Documentation: https://github.com/lh3/bwa

# This was the shell script I ran:

#!/bin/bash

#SBATCH --job-name=bwa
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load gcc/11.2.0
module load bwa/0.7.17
module load samtools/1.17

bwa index woodrat_hap1_purged.fa

bwa mem -t 128 woodrat_hap1_purged.fa combined_woodrat_L_R1.fastq.gz > hap1_purged_HiC_forward.sam
samtools sort -n -@128 -o hap1_purged_HiC_forward.bam hap1_purged_HiC_forward.sam

bwa mem -t 128 woodrat_hap1_purged.fa combined_woodrat_L_R2.fastq.gz > hap1_purged_HiC_reverse.sam
samtools sort -n -@128 -o hap1_purged_HiC_reverse.bam hap1_purged_HiC_reverse.sam

# This finished and took roughly an hour and a half. The resulting files were:

woodrat_hap1_purged.fa.bwt	
2.25 GB	8/5/2026 3:50:29 PM
file	woodrat_hap1_purged.fa.amb	
17 Bytes	8/5/2026 3:50:36 PM
file	woodrat_hap1_purged.fa.ann	
5.82 KB	8/5/2026 3:50:36 PM
file	woodrat_hap1_purged.fa.pac	
577 MB	8/5/2026 3:50:36 PM
file	woodrat_hap1_purged.fa.sa	
1.13 GB	8/5/2026 3:58:58 PM
file	hap1_purged_HiC_forward.sam	
350 GB	8/5/2026 5:12:09 PM
file	hap1_purged_HiC_forward.bam	
46.4 GB	8/5/2026 5:25:34 PM
file	hap1_purged_HiC_reverse.sam	
352 GB	8/5/2026 6:43:27 PM
file	slurm-19690523.out	
32.7 KB	8/5/2026 6:52:45 PM
file	hap1_purged_HiC_reverse.bam	
46.9 GB	8/5/2026 6:58:06 PM

# At this point the VGP does not give specifics to how they are filtering and merging and mapping. I refer here to YAHS to identify what the best
# inputs are as this is the scaffolder I will use


#---------------------------------------------------------------------------------------------------------------------------------------------

# Taking an alternate route here. I now plan to complete the arima hic mapping steps instead. This should accomplish the same steps as the 
# pre-processing step and is directly catered to our hic data.

#First I cloneed the arima HIC mapping pipeline:

git clone https://github.com/ArimaGenomics/mapping_pipeline.git

# The pipeline assumes the fasta is indexed so in the terminal I ran:

module load biocontainers/default
module load samtools/1.17

samtools faidx woodrat_hap1_purged.fa

#This produced the following, we will use this in the big script:

woodrat_hap1_purged.fa.fai

# Next I brought in the following files:

combined_woodrat_L_R1.fastq.gz	
55.4 GB	7/15/2026 12:44:52 PM

combined_woodrat_L_R2.fastq.gz	
55.6 GB	7/15/2026 12:46:09 PM

woodrat_hap1_purged.fa	
2.25 GB	8/11/2026 11:14:23 AM

# Then I modified the script from the pipeline for our file locations and the more recent tools:

#!/bin/bash

#####################################################
# ARIMA GENOMICS MAPPING PIPELINE Adapted 8/11/2026 #
#####################################################

#SBATCH --job-name=Arima_Mapping
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17
module load samtools/1.17
module load perl/5.32.0
module load picard/3.1.1

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw -p bwaindex/woodrat_hap1_purged woodrat_hap1_purged.fa

echo "### Step 1.A: FASTQ to BAM (R1)"
bwa mem -t 128 bwaindex/woodrat_hap1_purged combined_woodrat_L_R1.fastq.gz | samtools view -@ 128 -b - > woodrat_hap1_R1.bam

echo "### Step 1.B: FASTQ to BAM (R2)"
bwa mem -t 128 bwaindex/woodrat_hap1_purged combined_woodrat_L_R2.fastq.gz | samtools view -@ 128 -b - > woodrat_hap1_R2.bam

echo "### Step 2.A: Filter 5' end (R1)"
samtools view -h woodrat_hap1_R1.bam | perl filter_five_end.pl | samtools view -b - > woodrat_hap1_R1_filtered.bam

echo "### Step 2.B: Filter 5' end (R2)"
samtools view -h woodrat_hap1_R2.bam | perl filter_five_end.pl | samtools view -b - > woodrat_hap1_R2_filtered.bam

echo "### Step 3.A: Pair reads & mapping quality filter"
perl two_read_bam_combiner.pl woodrat_hap1_R1_filtered.bam woodrat_hap1_R2_filtered.bam samtools 10 | samtools view -b -t woodrat_hap1_purged.fa.fai - | samtools sort -@ 128 -o woodrat_hap1_HiC_paired.bam -

echo "### Step 3.B: Add read group"
picard AddOrReplaceReadGroups \
INPUT=woodrat_hap1_HiC_paired.bam \
OUTPUT=woodrat_hap1_HiC_RG.bam \
ID=woodrat_hap1_HiC \
LB=woodrat_hap1_HiC \
SM=woodrat \
PL=ILLUMINA \
PU=none

echo "### Step 4: Mark duplicates"
picard MarkDuplicates \
INPUT=woodrat_hap1_HiC_RG.bam \
OUTPUT=woodrat_hap1_HiC_dedup.bam \
METRICS_FILE=woodrat_hap1_HiC_dedup_metrics.txt \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT

samtools index woodrat_hap1_HiC_dedup.bam

perl get_stats.pl woodrat_hap1_HiC_dedup.bam > woodrat_hap1_HiC_dedup.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"


# This ran up to step 3A, perl is seemingly having trouble invoking samtools. To Fix this I created a wrapper 

#file name is samtools_wrapper.sh
# File contents:
#!/bin/bash
/usr/bin/singularity run /apps/biocontainers/images/alexanrna_samtools:1.17.sif samtools "$@"

# I then ran 
chmod +x samtools_wrapper.sh # This gives wrapper permission to run 

#To test that the wrapper works we will run a test 3A before rerunning:

perl two_read_bam_combiner.pl woodrat_hap1_R1_filtered.bam woodrat_hap1_R2_filtered.bam ./samtools_wrapper.sh 10 > combiner_test.sam

# This worked for the first step. I think I am going to rerun the entire arima mapping pipeline tonight with the fix:

#!/bin/bash

#####################################################
# ARIMA GENOMICS MAPPING PIPELINE Adapted 8/11/2026 #
#####################################################

#SBATCH --job-name=Arima_Mapping
#SBATCH --partition=highmem
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --time=24:00:00

cd $SLURM_SUBMIT_DIR

module load biocontainers/default
module load bwa/0.7.17
module load samtools/1.17
module load perl/5.32.0
module load picard/3.1.1

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw -p bwaindex/woodrat_hap1_purged woodrat_hap1_purged.fa

echo "### Step 1.A: FASTQ to BAM (R1)"
bwa mem -t 128 bwaindex/woodrat_hap1_purged combined_woodrat_L_R1.fastq.gz | samtools view -@ 128 -b - > woodrat_hap1_R1.bam

echo "### Step 1.B: FASTQ to BAM (R2)"
bwa mem -t 128 bwaindex/woodrat_hap1_purged combined_woodrat_L_R2.fastq.gz | samtools view -@ 128 -b - > woodrat_hap1_R2.bam

echo "### Step 2.A: Filter 5' end (R1)"
samtools view -h woodrat_hap1_R1.bam | perl filter_five_end.pl | samtools view -b - > woodrat_hap1_R1_filtered.bam

echo "### Step 2.B: Filter 5' end (R2)"
samtools view -h woodrat_hap1_R2.bam | perl filter_five_end.pl | samtools view -b - > woodrat_hap1_R2_filtered.bam

echo "### Step 3.A: Pair reads & mapping quality filter"
perl two_read_bam_combiner.pl woodrat_hap1_R1_filtered.bam woodrat_hap1_R2_filtered.bam ./samtools_wrapper.sh 10 | samtools view -b -t woodrat_hap1_purged.fa.fai - | samtools sort -@ 128 -o woodrat_hap1_HiC_paired.bam -
  
echo "### Step 3.B: Add read group"
picard AddOrReplaceReadGroups \
INPUT=woodrat_hap1_HiC_paired.bam \
OUTPUT=woodrat_hap1_HiC_RG.bam \
ID=woodrat_hap1_HiC \
LB=woodrat_hap1_HiC \
SM=woodrat \
PL=ILLUMINA \
PU=none

echo "### Step 4: Mark duplicates"
picard MarkDuplicates \
INPUT=woodrat_hap1_HiC_RG.bam \
OUTPUT=woodrat_hap1_HiC_dedup.bam \
METRICS_FILE=woodrat_hap1_HiC_dedup_metrics.txt \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT

samtools index woodrat_hap1_HiC_dedup.bam

perl get_stats.pl woodrat_hap1_HiC_dedup.bam > woodrat_hap1_HiC_dedup.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"



















