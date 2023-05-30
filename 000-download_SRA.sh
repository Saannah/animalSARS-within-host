Very useful link:
https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump

Nodes on Narval do not have internet access and only the login node has internet. I used prefetch to download the files into Narval. 

1. load modules on login node
2. open configuration to change directory and disconnect external network
3. download using prefetch
4. validate using vdb-validate
5. convert to fastq paired reads using fastqdump and --split-3 option

########### script to download fastq after prefetching
#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a job submission on a Compute Canada cluster.
# ---------------------------------------------------------------------
#SBATCH --account=def-shapiro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=10:00  # e.g. 1-00:00:00 for a day
#SBATCH --ntasks=1
#SBATCH --output=dl.log


module load StdEnv/2020  gcc/9.3.0
module load sra-toolkit/3.0.0
cd $1
fastq-dump --split-3 ${1}.sra