#!/bin/bash -l

#SBATCH --job-name=M-pre-qc
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/.../merge/MD/log/7a_pre_qc_imp.out
#SBATCH --error=/gpfs/.../merge/MD/log/7a_pre_qc_imp.err

# load modules
module load gcc
module load plink/1.9
module load vcftools
module load R
module load bcftools

#Set arguments
plink_prefix=MDs_forImp

out_dir=/gpfs/.../merge/MD

cd  ${out_dir}
bash 7b_create_initial_input.sh ${plink_prefix} ${out_dir}
