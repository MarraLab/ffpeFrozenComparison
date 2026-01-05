#!/bin/bash -l

#SBATCH --job-name=batchcorrect # Job name
#SBATCH --nodes=1 # Nodes requested
#SBATCH --cores=10 # Cores requested
#SBATCH --mem=200G # max memory allocated
#SBATCH --output=/projects/marralab/scratch/cayan/techComp/%x.o.%j
#SBATCH --mail-type=ALL # email if failed
#SBATCH --mail-user=cayan # email me
#SBATCH --error=/projects/marralab/scratch/cayan/techComp/%x.e.%j

# source /home/cayan/.bashrc
source ~/anaconda3/etc/profile.d/conda.sh
# conda activate scenv

# Ensure Conda libraries are visible to R
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/cayan/anaconda3/lib

# input=$1

/gsc/software/linux-x86_64-centos7/R-4.2.2/bin/Rscript batchCorrectEval_function.R "$input"