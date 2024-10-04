#!/bin/bash -l

#SBATCH --job-name=M-pre-qc
#SBATCH --partition=tier1q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/data/.../merge/MD/log/8fix_strands.out
#SBATCH --error=/gpfs/data/.../merge/MD/log/8fix_strands.err

# load modules
module load gcc/12.1.0
module load plink/1.9
module load vcftools
module load R/4.3.1
module load bcftools/1.17


# since bcftools +fixref fliped strands so that this step could be skipped 
# no snps excluded due to 'Strand flip' so 4get_strand_flip_snp_names.R doesn't work

#Set arguments
pre_qc_dir=/gpfs/.../merge/MD/pre_qc
post_qc_dir=/gpfs/.../merge/MD/post_qc

bash 8a_fix_strands.sh ${pre_qc_dir} ${post_qc_dir}
