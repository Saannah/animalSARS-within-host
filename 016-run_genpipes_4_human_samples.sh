#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a job submission on a Compute Canada cluster.
# ---------------------------------------------------------------------
#SBATCH --account=ctb-shapiro
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1-00:00:00  # e.g. 1-00:00:00 for a day
#SBATCH --ntasks=1
#SBATCH --output=genpipes


export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/python/3.9.1
module load mugqic/genpipes/4.3.0
export JOB_MAIL=sana.naderi@mail.mcgill.ca
export RAP_ID=ctb-shapiro

covseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/covseq/covseq.base.ini \
$MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini \
$MUGQIC_INSTALL_HOME/genomes/species/Coronavirinae.SARS-CoV-2/Coronavirinae.SARS-CoV-2.ini \
-j slurm \
-r human_readset.txt \
-g out.sh

$MUGQIC_PIPELINES_HOME/utils/chunk_genpipes.sh -n 20 out.sh chunks

$MUGQIC_PIPELINES_HOME/utils/submit_genpipes chunks
