#!/bin/bash
#
#SBATCH --job-name=compute_prevalence
#SBATCH --output=/scratch/users/cdonnat/aerosol_transmission_model/out/prevalence.out
#SBATCH --error=/scratch/users/cdonnat/aerosol_transmission_model/out/error.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --qos=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners
# load modules

ml R
for ind in {1..192}
do
   sbatch bash_script ${ind}
done
