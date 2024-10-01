#!/bin/bash -l

#SBATCH --job-name=subset
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=tier2q
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/.../log/5subset.out
#SBATCH --error=/gpfs/.../log/5subset.err

module load gcc
module load plink

plink_prefix=MD2123_het11
out_dir=/gpfs/.../MD2123
cd ${out_dir}

# create samples list
bash extrSmpByStudy.sh

# Convert to vcf

plink2 \
        --bfile ${out_dir}/${plink_prefix} \
        --keep subset.list \
        --export vcf id-paste=iid \
        --threads 14 \
        --out ${out_dir}/subset/subset_QCd


