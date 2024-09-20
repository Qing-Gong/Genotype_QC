#!/bin/bash -l

#SBATCH --job-name=GDA_convert
#SBATCH --nodes=4
#SBATCH --partition=tier3q
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/data/scratch/users/qgong/MP/log/1MD.out
#SBATCH --error=/gpfs/data/scratch/users/qgong/MP/log/1MD.err

module load gcc
module load plink/1.9

# set working directory 
cd /gpfs/data/scratch/users/qgong/MP

# convert Illumina report to plink files
SAMPLE_MAP=/gpfs/data/lab/resources/MP/genotype_data/Illumina_GDAv8/Excel_Files/MP_Sample_Map.txt
SAMPLE_MAP_FILE=$(basename ${SAMPLE_MAP})
PREFIX=$(echo ${SAMPLE_MAP_FILE}|sed 's/_Sample_Map\.txt//')

REL_PATH=$(dirname ${SAMPLE_MAP})

PLATE_PATH=/gpfs/data/lab/resources/MP/genotype_data

awk 'NR < 2 { next }  !seen[$2]++ {print $2,$2,0,0,$3,0}' ${SAMPLE_MAP} > ${PREFIX}.fam 

awk 'NR < 2 { next } {print $3,$2,0,$4}' ${REL_PATH}/${PREFIX}_SNP_Map.txt > ${PREFIX}.map

awk 'NR < 11 { next }  {print $2,$2,$1,$3,$4}' ${REL_PATH}/${PREFIX}_FinalReport.txt > ${PREFIX}.lgen


# REFORMAT
sed 's/- -/0 0/g' ${PREFIX}.lgen > ${PREFIX}_fmt.lgen
sed 's/Male/1/g' ${PREFIX}.fam > ${PREFIX}_fmt.fam
sed -i 's/Female/2/g' ${PREFIX}_fmt.fam
sed -i 's/Unknown/0/g' ${PREFIX}_fmt.fam

mv ${PREFIX}.map ${PREFIX}_fmt.map

## CONVERT TO PLINK:

plink --lfile ${PREFIX}_fmt --recode --out ${PREFIX}_plink
plink --lfile ${PREFIX}_fmt --make-bed --out ${PREFIX}_plink

# Update sampleI D
#sed '1,22d' ${PLATE_PATH}/MP_2_genotyping_plate_map.csv|awk -F"," '{printf("%s\t%s_%s\n",$2, $1, "p15")}' >plate15.txt

sed '1d;25,31d;33d' ${PLATE_PATH}/MP_plate_map.csv|awk -F"," '{printf("%s\t%s_%s\n",$2, $1, "p15")}' >plate15.txt

paste <(awk '{print $1, $2}' ${PREFIX}_fmt.fam) <(awk '{print $2,$2}' plate15.txt) >p15_sample.txt

plink \
        --bfile ${PREFIX}_plink \
        --keep-allele-order \
        --update-ids p15_sample.txt \
        --make-bed \
        --out MP-p15

# Extract SUBSET
plink \
	--bfile MP-p15 \
	--remove MP21.list \
	--recode \
	--out HG_p15

# Strand correction
# step 1: illumina manifest file checking
/gpfs/data/ -lab/resources/TOOLS/Strandscript/bin/step1-mismatch.pl \
        -g hg19 \
        -in /gpfs/data/ -scratch/users/qgong/reference/illumina/GDA-8v1-0_D1.csv \
        -n GDA

# step 2: strand flip for plink file (input:ped & map)
/gpfs/data/ -lab/resources/TOOLS/Strandscript/bin/step2-flip.pl \
        -in new_GDA.csv \
        -map HG_p15.map \
        -ped HG_p15.ped

# remove duplicates
plink \
        --file flipped_HG_p15 \
        --list-duplicate-vars ids-only suppress-first \
        --out plink

plink \
	--file flipped_HG_p15 \
	--exclude plink.dupvar \
        --snps-only just-acgt \
        --make-bed \
        --out HG_0dup0

# Update sex and remove SNPs on chromosome 0
plink \
        --bfile HG_0dup0 \
        --update-sex HG_p15.sex \
        --not-chr 0 \
        --make-bed \
        --out HG_updSex1

# check Sex and stats
plink \
    --bfile HG_updSex1 \
    --check-sex \
    --missing \
    --freq \
    --hardy \
    --out plink_stat
