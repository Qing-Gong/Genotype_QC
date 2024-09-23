#!/bin/bash -l

#SBATCH --job-name=liftOver
#SBATCH --nodes=1
#SBATCH --partition=tier1q
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/data/scratch/users/qgong/MD/log/3liftOver.out
#SBATCH --error=/gpfs/data/scratch/users/qgong/MD/log/3liftOver.err

module load gcc
module load plink/1.9
module load python/3.10.5
module load R


# Create sex file based on R script
Rscript 3a_MD_sex_race.R

sed '1d' MD_sex_race_AFR.csv |awk -F',' '{printf "%s\t%s\t%s\n", $1,$1,$2}' >MD_AMRAFR.sex

# Update sex
plink \
    --bfile MD2123_0dup0 \
    --update-sex MD_AMRAFR.sex \
    --not-chr 0 \
    --make-bed \
    --out MD2123_updSex1

# check Sex and stats
plink \
    --bfile MD2123_updSex1 \
    --check-sex \
    --missing \
    --freq \
    --hardy \
    --out plink_stat

# plink2 changed 23, 24, 25 and 26 to X,Y, XY and MT

# Create bed file to crossover from hg19 to hg38 
awk '{print "chr" $1, $4, $4, $2}' MD2123_updSex1.bim|sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > MD2123_updSex1.tolift

# Do crossover
CrossMap.py bed /gpfs/data/scratch/users/qgong/reference/liftOver/hg19ToHg38.over.chain \
      MD2123_updSex1.tolift  \
      MD2123_updSex1_hg38

# Extract mapped variants
awk '{print $4}' MD2123_updSex1_hg38 > MD2123_liftOver.snps

plink \
    --bfile MD2123_updSex1 \
    --extract MD2123_liftOver.snps \
    --make-bed \
    --out MD2123_lift2

# Update positions   # avoid using --update-map and getting warnig "Base-pair positions are now unsorted!"
cat /gpfs/data/scratch/users/qgong/MD/3b_update_pos.R | R --vanilla

# king -b MD2123_lift2.bed --kinship



