**On the login node, load the modules:
module load StdEnv/2020  gcc/9.3.0
module load sra-toolkit/3.0.0



https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
**run the configuration using the following command and do the following:
vdb-config -i
1. in tools change working directory
2. in main, make sure make sure that Enable Remote Access is not checked.


**download all the accessions usinng prefetch:
for srr in $(cat accession_mink.txt)
do
prefetch -O /home/sanna/projects/def-shapiro/sanna/aim3/sra/mink $srr
done



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






**** call the above script in a loop

