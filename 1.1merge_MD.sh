#!/bin/bash

#PBS -N merge
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=04:00:00
#PBS -e /scratch/qgong/MD/log/merge.err
#PBS -o /scratch/qgong/MD/log/merge.out

## CONVERT TO PLINK:
module load gcc/6.2.0
module load plink/1.90

cd /scratch/qgong/MOD-2022-08/mod2021/getPlink

plink \
	--bfile /scratch/qgong/MD/MD2023 \
	--merge-list /scratch/qgong/MD/scripts/MD2123/mergeFiles_update.txt \
	--keep-allele-order \
	--check-sex \
	--missing \
	--freq \
	--hardy \
	--list-duplicate-vars ids-only suppress-first \
	--make-bed \
	--out MD2123
