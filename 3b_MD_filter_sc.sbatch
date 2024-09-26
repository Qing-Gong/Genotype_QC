#!/bin/bash -l

#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=32gb
#SBATCH --error=/scratch/qgong/.../log/3b_MD_flt.err
#SBATCH --output=/scratch/qgong/.../log/3b_MD_flt.out

#Set arguments
plink_prefix=MD2123
out_dir=/scratch/qgong/.../getPlink
script_dir=/scratch/qgong/.../scripts/MD2123

cd ${out_dir}

MDule load gcc/11.3.0
MDule load plink/1.9
MDule load R
MDule load king

# remove monomorphic SNPs, SNPs with high
# missingness, SNPs not in HWE
plink \
        --bfile MD2123_liftOver_no_AT_CG3 \
        --geno 0.03 \
	--check-sex \
        --make-bed \
        --out ${plink_prefix}_geno4

plink \
        --bfile ${plink_prefix}_geno4 \
        --mind 0.03 \
        --missing \
        --make-bed \
        --out ${out_dir}/MD2123_mind5

plink --bfile ${plink_prefix}_mind5 \
      --remove ${script_dir}/non_MD.csv \
      --make-bed \
      --out ${plink_prefix}_rm_nonMD6

king -b MD2123_rm_nonMD6.bed --kinship

# Find replicates and samples need correcting ID
# Thanks WM
python3 ${script_dir}/gen_report.py -m ${plink_prefix}_mind5.imiss -k king.kin0 -c ${script_dir}/MD_samples_change_2123.txt

# Create samples need removing
less MD_ids_remove.txt|awk '{printf "%s\t%s\n", $1, $1}' > bad_samples.txt

plink \
	--bfile ${plink_prefix}_rm_nonMD6 \
	--remove bad_samples.txt \
	--make-bed \
	--out ${plink_prefix}_rm_badID7
 
plink \
	--bfile ${plink_prefix}_rm_badID7 \
	--update-ids MD_ids_relabel.txt \
        --make-bed \
        --out ${plink_prefix}_crtID8
less ${plink_prefix}_crtID8.fam|awk '{print $1}' >MD_crtID8.list

# Seperate samples by race and by mom-babies
Rscript --no-save ${script_dir}/3b_separate_race_kin.R

# Filter by HWE by race and mom/babies
plink \
        --bfile ${plink_prefix}_crtID8 \
        --recode \
        --keep black_baby.txt \
        --hardy midp \
        --out qc_hwe_blackB

plink \
        --bfile ${plink_prefix}_crtID8 \
        --recode \
        --keep black_mom.txt \
        --hardy midp \
        --out qc_hwe_blackM

plink \
        --bfile ${plink_prefix}_crtID8 \
        --recode \
        --keep white_baby.txt \
        --hardy midp \
        --out qc_hwe_whiteB

plink \
        --bfile ${plink_prefix}_crtID8 \
        --recode \
        --keep white_mom.txt \
        --hardy midp \
        --out qc_hwe_whiteM

# get SNPs violating HWE according to meta-p < 5*10-8 in each reported race and mom/baby

awk '{if ($9<1e-3) {print $2}}' qc_hwe_blackB.hwe > bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_blackM.hwe >> bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_whiteB.hwe >> bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_whiteM.hwe >> bad_HWE_snps.list

# filter out bad HWE SNPs
plink \
        --bfile ${plink_prefix}_crtID8 \
        --exclude bad_HWE_snps.list \
        --het \
        --make-bed  \
        --out ${plink_prefix}_hwe9

Rscript --no-save ${script_dir}/heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# remove samples failed het
plink \
        --bfile ${plink_prefix}_hwe9 \
        --remove het_fail_ind.txt \
        --make-bed \
        --out ${plink_prefix}_het10
