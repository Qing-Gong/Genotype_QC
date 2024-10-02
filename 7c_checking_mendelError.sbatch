#!/bin/bash -l

#SBATCH --job-name=Mendel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/.../merge/log/3checkingME.out
#SBATCH --error=/gpfs/.../merge/log/3checkingME.err

module load gcc
module load plink/1.9
module load R
module load king

plink_prefix=MRNH
out_dir=/gpfs/.../merge
cd ${out_dir}

Rscript 3a_mom_baby_pair.R

<MD_mom_baby_pair.csv tr -s ' ' '\n' >mom_baby.list

awk '{printf "%s\t%s\n", 0, $1}' mom_baby.list >moms_babies.txt 

# check stat
plink \
        --bfile ${plink_prefix}_hwe3 \
        --threads 14 \
	--keep moms_babies.txt \
	--freq \
	--make-bed \
	--out ${plink_prefix}_mb4

king -b ${plink_prefix}_mb4.bed --kinship --prefix king3

Rscript 3b_addMomID.R

cp ${plink_prefix}_mb4.bim ${plink_prefix}_mb4_update.bim
cp ${plink_prefix}_mb4.bed ${plink_prefix}_mb4_update.bed

plink \
	--bfile ${plink_prefix}_mb4_update \
	--mendel \
	--mendel-duos \
        --out Medelian_${plink_prefix}

plink \
	--bfile ${plink_prefix}_mb4_update \
        --me 0.05 0.1 \
	--mendel-duos \
        --make-bed \
        --out ${plink_prefix}_mendel5

#        --mendel-duos \
#	--set-me-missing \
#	--make-bed \
#	--out ${plink_prefix}_mendel5
