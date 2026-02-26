#!/bin/bash -l

#SBATCH --job-name=batchcorrect # Job name
#SBATCH --nodes=1 # Nodes requested
#SBATCH --cores=10 # Cores requested
#SBATCH --mem=200G # max memory allocated
#SBATCH --output=output
#SBATCH --mail-type=ALL # email if failed
#SBATCH --mail-user=cayan # email me
#SBATCH --error=error

# source /home/cayan/.bashrc
source ~/anaconda3/etc/profile.d/conda.sh
# conda activate scenv

# Ensure Conda libraries are visible to R
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path

# input=$1

Rscript batchCorrectEval_function.R "$input"