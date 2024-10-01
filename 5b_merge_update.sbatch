#!/bin/bash -l

#SBATCH --job-name=merge1
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --partition=tier2q
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=64gb
#SBATCH --output=/gpfs/.../merge/log/1merge_update.out
#SBATCH --error=/gpfs/.../merge/log/1merge_update.err

module load gcc
module load htslib
module load bcftools
module load plink

out_dir=/gpfs/.../merge
cd ${out_dir}


# AMR-AFR
plink_prefix1=MD2123
dir1=/gpfs/.../MD2123

# Convert to vcf
plink2 \
        --bfile ${dir1}/${plink_prefix1}_het11 \
        --export vcf id-paste=iid \
        --threads 14 \
        --out ${plink_prefix1}

# rename 1 to chr1
<${plink_prefix1}.vcf sed 's/^\([0-9]\)/chr\1/;s/^X/chrX/;s/^Y/chrY/' >${plink_prefix1}_chr.vcf

# MEGA
plink_prefix2=MD2022
dir2=/gpfs/.../MD2022

# Convert to vcf
plink2 \
        --bfile ${dir2}/${plink_prefix2}_het11 \
        --threads 14 \
        --export vcf id-paste=iid \
        --out ${plink_prefix2}

# rename 1 to chr1
<${plink_prefix2}.vcf sed 's/^\([0-9]\)/chr\1/;s/^X/chrX/;s/^Y/chrY/' >${plink_prefix2}_chr.vcf

# GDA
plink_prefix3=HM
dir3=/gpfs/.../HG2023

# Convert to vcf
plink2 \
        --bfile ${dir3}/${plink_prefix3}_het8 \
        --threads 14 \
        --export vcf id-paste=iid \
        --out ${plink_prefix3}

<${plink_prefix3}.vcf sed 's/^\([0-9]\)/chr\1/;s/^X/chrX/;s/^Y/chrY/' >${plink_prefix3}_chr.vcf


# Compress the vcf file
for files in ${plink_prefix1}_chr.vcf ${plink_prefix2}_chr.vcf ${plink_prefix3}_chr.vcf
        do
                bgzip "$files"
        done

# Build index and add annotation
for files in ${plink_prefix1}_chr.vcf.gz ${plink_prefix2}_chr.vcf.gz ${plink_prefix3}_chr.vcf.gz
        do
                tabix -p vcf "$files"
        done

# consistent reference SNPs--flip/swap the reference and alternate alleles  -m flip not flip ambiguous sites but -m flip-all  for flipping any sites, including ambiguous A/T and C/G pairs
bcftools +fixref ${plink_prefix1}_chr.vcf.gz -- -f /gpfs/data/referenceFiles/Homo_sapiens/STARgenome/GRCh38.primary_Gencode38/Sequence/GRCh38.primary_assembly.genome.fa -i /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -m flip-all >${plink_prefix1}_ref.vcf
bcftools +fixref ${plink_prefix2}_chr.vcf.gz -- -f /gpfs/data/referenceFiles/Homo_sapiens/STARgenome/GRCh38.primary_Gencode38/Sequence/GRCh38.primary_assembly.genome.fa -i /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -m flip-all >${plink_prefix2}_ref.vcf
bcftools +fixref ${plink_prefix3}_chr.vcf.gz -- -f /gpfs/data/referenceFiles/Homo_sapiens/STARgenome/GRCh38.primary_Gencode38/Sequence/GRCh38.primary_assembly.genome.fa -i /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -m flip-all >${plink_prefix3}_ref.vcf


# The above command might have changed the coordinates, we must sort the VCF.
bcftools sort ${plink_prefix1}_ref.vcf -Ov -o ${plink_prefix1}_ref_sorted.vcf
bcftools sort ${plink_prefix2}_ref.vcf -Ov -o ${plink_prefix2}_ref_sorted.vcf
bcftools sort ${plink_prefix3}_ref.vcf -Ov -o ${plink_prefix3}_ref_sorted.vcf


# Compress the vcf file
for files in ${plink_prefix1}_ref_sorted.vcf ${plink_prefix2}_ref_sorted.vcf ${plink_prefix3}_ref_sorted.vcf
        do
                bgzip "$files"
        done

# Build index and add annotation
for files in ${plink_prefix1}_ref_sorted.vcf.gz ${plink_prefix2}_ref_sorted.vcf.gz ${plink_prefix3}_ref_sorted.vcf.gz
        do
                tabix -p vcf "$files"
        done


# add annotation 
bcftools annotate -a /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -c ID -o ${plink_prefix1}_anno.vcf ${plink_prefix1}_ref_sorted.vcf.gz

bcftools annotate -a /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -c ID -o ${plink_prefix2}_anno.vcf ${plink_prefix2}_ref_sorted.vcf.gz

bcftools annotate -a /gpfs/.../ref/dbsnp/dbSNP155.hg38.gz -c ID -o ${plink_prefix3}_anno.vcf ${plink_prefix3}_ref_sorted.vcf.gz

# Compress the vcf file
for files in ${plink_prefix1}_anno.vcf ${plink_prefix2}_anno.vcf ${plink_prefix3}_anno.vcf
        do
                bgzip "$files"
        done

# Build index and add annotation
for files in ${plink_prefix1}_anno.vcf.gz ${plink_prefix2}_anno.vcf.gz ${plink_prefix3}_anno.vcf.gz
        do
                tabix -p vcf "$files"
        done


# Determine SNPs in common
# Determine SNPs in common -c none only records with identical REF and ALT alleles are compatible -all all records are compatible, regardless of whether the ALT alle$

bcftools isec -p ${out_dir} -n=3 -c none ${plink_prefix1}_anno.vcf.gz ${plink_prefix2}_anno.vcf.gz ${plink_prefix3}_anno.vcf.gz

#compress vcf files
bgzip 0000.vcf
bgzip 0001.vcf
bgzip 0002.vcf

#build index
tabix -p vcf 0000.vcf.gz
tabix -p vcf 0001.vcf.gz
tabix -p vcf 0002.vcf.gz

##-m none .. no new multiallelics, output multiple records instead -id merge by ID
bcftools merge -m id 0000.vcf.gz 0001.vcf.gz 0002.vcf.gz -o MRNH.vcf


