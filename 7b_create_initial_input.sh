#!/bin/bash

#Set arguments
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <plink_file_prefix> <out_dir>"
    exit
fi
plink_prefix=$1
out_dir=$2

#Remove SNPs with duplicate positions
plink --bfile ${plink_prefix} \
   --list-duplicate-vars suppress-first \
   --out tmp_dupl_check

cat tmp_dupl_check.dupvar | sed -e '1d' | cut -f4 > tmp_dupl_snpids.txt

plink --bfile ${plink_prefix} \
   --exclude tmp_dupl_snpids.txt \
   --make-bed --out tmp_no_dupl

#Remove strand ambiguous SNPs
cat 3get_strand_amb_SNPs.R | R --vanilla

plink --bfile tmp_no_dupl \
      --exclude tmp_strand_remove_snps.txt \
      --make-bed --out tmp_no_AT_CG

#Perform pre-imputation QC - remove monomorphic SNPs, SNPs with high
#missingness, SNPs not in HWE  # --maf 0.000001 

plink \
        --bfile tmp_no_AT_CG \
        --maf 0.000001 \
        --geno 0.05 \
        --hwe 5e-8 \
        --make-bed \
        --out ${out_dir}/pre_qc/pre_qc

#Create vcf files for uploading to imputation server for QC
#Note that the encoding for chromosome is e.g. chr22, not chr
for ((chr=1; chr<=22; chr++)); do
    plink --bfile ${out_dir}/pre_qc/pre_qc --chr $chr --recode vcf --out tmp_chr${chr}
    vcf-sort tmp_chr${chr}.vcf | sed -E 's/^([[:digit:]]+)/chr\1/' | bgzip -c > ${out_dir}/pre_qc/chr${chr}_pre_qc.vcf.gz
done

#Report SNP counts
orig_snp_nr=`wc -l ${plink_prefix}.bim`
nonamb_snp_nr=`wc -l tmp_no_AT_CG.bim`
qc_snp_nr=`wc -l ${out_dir}/pre_qc/pre_qc.bim`
echo "Original SNP nr: $orig_snp_nr"
echo "Non-ambiguous SNP nr: $nonamb_snp_nr"
echo "Final SNP nr after QC: $qc_snp_nr"

#Cleanup
rm tmp_*
