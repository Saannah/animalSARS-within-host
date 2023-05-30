### how to run Freyja using Steven's singularity container:


## load singularity
module load singularity/3.8

## store all bam files + the .sif container + reference fasta in a directory (004-freyja)
## go to the directory above (../) and run the following:

singularity exec --no-home -B 004-freyja/:/mnt 004-freyja/Freyja_V6.sif freyja variants /mnt/SRR14095171.sorted.filtered.bam --variants /mnt/variant_out.tsv --depths /mnt/depth_file.tsv --ref /mnt/sars-cov2-ref-used-by-genpipes.fasta
singularity exec --no-home -B 004-freyja/:/mnt 004-freyja/Freyja_V6.sif freyja demix /mnt/variant_out.tsv /mnt/depth_file.tsv --output /mnt/output.tsv



#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a job submission on a Compute Canada cluster.
# ---------------------------------------------------------------------
#SBATCH --account=ctb-shapiro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00  # e.g. 1-00:00:00 for a day
#SBATCH --ntasks=1
#SBATCH --output=log.freyja

module load singularity/3.8
for i in *bam
do
   singularity exec --no-home -B run_w_singularity/:/mnt run_w_singularity/Freyja_V6.sif freyja variants /mnt/${i} --variants /mnt/${i}_variant.tsv --depths /mnt/${i}_depth_file.tsv --ref /mnt/sars-cov2-ref-used-by-genpipes.fasta
   singularity exec --no-home -B run_w_singularity/:/mnt run_w_singularity/Freyja_V6.sif freyja demix /mnt/${i}_variant.tsv /mnt/${i}_depth_file.tsv --output /mnt/${i}_output.tsv
done
