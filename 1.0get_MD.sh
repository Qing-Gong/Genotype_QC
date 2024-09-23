#!/bin/bash
#PBS -N extract5 MD_AFR
#PBS -l walltime=00:04:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=32gb 
#PBS -e /scratch/qgong/MD/log/MD_AFR5.err
#PBS -o /scratch/qgong/MD/log/MD_AFR5.out

# set working directory 
cd /scratch/qgong/MD/getPlink

# convert Illumina report to plink files
SAMPLE_MAP=/gpfs/data/lab/resources/.../Excel_Data/MD_Sample-Map.txt
SAMPLE_MAP_FILE=$(basename ${SAMPLE_MAP})
PREFIX=$(echo ${SAMPLE_MAP_FILE}|sed 's/_Sample-Map\.txt//')

REL_PATH=$(dirname ${SAMPLE_MAP})

PLATE_PATH=/scratch/qgong/MD/plates_submission


awk 'NR < 2 { next }  !seen[$2]++ {print $2,$2,0,0,$3,0}' ${SAMPLE_MAP} > ${PREFIX}.fam 

awk 'NR < 2 { next } {print $3,$2,0,$4}' ${REL_PATH}/${PREFIX}_SNP-Map.txt > ${PREFIX}.map

awk 'NR < 11 { next }  {print $2,$2,$1,$3,$4}' ${REL_PATH}/${PREFIX}_FinalReport.txt > ${PREFIX}.lgen


# REFORMAT
sed 's/- -/0 0/g' ${PREFIX}.lgen > ${PREFIX}_fmt.lgen
sed 's/Male/1/g' ${PREFIX}.fam > ${PREFIX}_fmt.fam
sed -i 's/Female/2/g' ${PREFIX}_fmt.fam
sed -i 's/Unknown/0/g' ${PREFIX}_fmt.fam

mv ${PREFIX}.map ${PREFIX}_fmt.map

MDule load gcc/6.2.0
MDule load plink/1.90

plink \
        --lfile ${PREFIX}_fmt \
        --recode \
        --out ${PREFIX}_plink

plink \
        --file ${PREFIX}_plink \
        --make-bed \
        --out ${PREFIX}_plink

#update sampleID
dos2unix ${PLATE_PATH}/Plate5_Plate_Map.csv      ##remove line-ending style charater(\r) CRLF
sed '1,1d' ${PLATE_PATH}/Plate5_Plate_Map.csv|head -96|awk -F"," '{gsub(/ /,"",$2);printf("%s\t%s_%s\n", $1, $2, "p5")}'>plate5.txt
paste <(awk '{print $1, $2}' ${PREFIX}_fmt.fam) <(awk '{print $2,$2}' plate5.txt) >p5_sample.txt

plink \
        --file ${PREFIX}_plink \
        --keep-allele-order \
        --update-ids p5_sample.txt \
        --list-duplicate-vars ids-only suppress-first \
        --make-bed \
        --out MD-p5

