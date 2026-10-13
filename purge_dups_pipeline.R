#Purge_dups pipeline

#In order to figure out how to run purge_dups I creeated this directory to trouble shoot how to gewt it to run correctly. These steps and text is adapted from the purge_dups pipeline on github
#https://github.com/dfguan/purge_dups#usg


#Pipeline

#Given a primary assembly pri_asm and an alternative assembly hap_asm (optional, if you have one), follow the steps shown below to build your own purge_dups pipeline, 
#steps with same number can be run simultaneously. Among all the steps, although step 4 is optional, we highly recommend our users to do so, because assemblers may produce 
#overrepresented seqeuences. In such a case, The final step 4 can be applied to remove those seqeuences.

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Step 1A. Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth. Commands are as follows:
#For PacBio CCS reads. Each Step will be run in its own directory seperated by haplotype subdirectories to prevent overwritting the same file names

#Step 1A. Hap_1

#!/bin/bash

#SBATCH --job-name=step1A_hap1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

minimap2 -xmap-hifi -I 4G woodrat_hap1.fasta m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz | gzip -c - > hap_1.paf.gzI 
pbcstat hap_1.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log


#Step 1A. Hap_2

#!/bin/bash

#SBATCH --job-name=step1A_hap2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

minimap2 -xmap-hifi -I 4G woodrat_hap2.fasta m84082_260311_175030_s3.hifi_reads.bc2024.fastq.gz | gzip -c - > hap_2.paf.gzI 
pbcstat hap_2.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

#Notice If you have a large genome, please set minimap2 -I option to ensure the genome can be indexed once, otherwise read depth can be wrong.

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
#Step 1B. Split an assembly and do a self-self alignment. Commands are following:

#Step 1B. Hap_1

#!/bin/bash

#SBATCH --job-name=step1B_hap1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

split_fa woodrat_hap1.fasta > woodrat_hap1.fasta.split
minimap2 -xasm5 -DP woodrat_hap1.fasta.split woodrat_hap1.fasta.split | gzip -c - > woodrat_hap1.fasta.split.self.paf.gz


#Step 1B. Hap_2

#!/bin/bash

#SBATCH --job-name=step1B_hap2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

split_fa woodrat_hap2.fasta > woodrat_hap2.fasta.split
minimap2 -xasm5 -DP woodrat_hap2.fasta.split woodrat_hap2.fasta.split | gzip -c - > woodrat_hap2.fasta.split.self.paf.gz

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Step 2. Purge haplotigs and overlaps with the following command.

#Step 2. Hap_1

#!/bin/bash

#SBATCH --job-name=step2_hap1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

purge_dups -2 -T cutoffs -c PB.base.cov woodrat_hap1.fasta.split.self.paf.gz > hap1_dups.bed 2> hap1_purge_dups.log


#Step 2. Hap_2

#!/bin/bash

#SBATCH --job-name=step2_hap2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

purge_dups -2 -T cutoffs -c PB.base.cov woodrat_hap2.fasta.split.self.paf.gz > hap2_dups.bed 2> hap2_purge_dups.log

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Step 3. Get purged primary and haplotig sequences from draft assembly.

#Step 3. Hap_1

#!/bin/bash

#SBATCH --job-name=step3_hap1
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

get_seqs -e hap1_dups.bed woodrat_hap1.fasta


#Step 3. Hap_2

#!/bin/bash

#SBATCH --job-name=step3_hap2
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00

module load biocontainers/default
module load purge_dups/1.2.6_fixed

get_seqs -e hap2_dups.bed woodrat_hap2.fasta

#Notice this command will only remove haplotypic duplications at the ends of the contigs. 
#If you also want to remove the duplications in the middle, please remove -e option at your own risk, 
#it may delete false positive duplications. For more options, please refer to get_seqs -h.

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Step 4. Merge hap.fa and $hap_asm and redo the above steps to get a decent haplotig set.


#Limitation
Read depth cutoffs calculation: the coverage cutoffs can be larger for a low heterozygosity species, which causes the purged assembly size smaller than expected. In such a case, please use script/hist_plot.py to make the histogram plot and set coverage cutoffs manually.
Repeats: purge_dups has a limited ability to process repeats.
  
  
  
  
  
  
  
  
  
  
  
  
  
  