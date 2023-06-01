#!/bin/bash

#SBATCH --account=smed003061
#SBATCH --job-name QC_GWAS
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem=80G
#SBATCH --time=08:00:00

cd /user/work/hd23261

export R_LIBS=Rlibs

module load languages/r/4.2.1

Rscript GWASinspector.R

