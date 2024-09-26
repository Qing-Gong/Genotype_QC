#!/bin/bash -l

#SBATCH --job-name=filter
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=32gb
#SBATCH --error=/gpfs/data/.../qgong/.../log/3MD2123_flt.err
#SBATCH --output=/gpfs/data/.../qgong/.../log/3MD2123_flt.out

#Set arguments
plink_prefix=MD2123
out_dir=/gpfs/data/.../qgong/.../MD2123

cd ${out_dir}

module load gcc
module load plink/1.9
module load R
module load king

# remove monomorphic SNPs, SNPs with high
# missingness, SNPs not in HWE
plink \
        --bfile ${plink_prefix}_lift2 \
        --geno 0.1 \
        --make-bed \
        --out ${plink_prefix}_geno3

plink \
        --bfile ${plink_prefix}_geno3 \
        --mind 0.1 \
        --make-bed \
        --out ${plink_prefix}_mind4

plink \
        --bfile ${plink_prefix}_mind4 \
        --geno 0.05 \
        --make-bed \
        --out ${plink_prefix}_geno5

plink \
        --bfile ${plink_prefix}_geno5 \
        --mind 0.05 \
        --missing \
        --check-sex \
        --make-bed \
        --out ${plink_prefix}_mind6

king -b ${plink_prefix}_mind6.bed --kinship 

# Find replicates and samples need correcting ID # Compared replicates pair and remove replicates with larger missingness rate

python3 gen_report.py -m ${plink_prefix}_mind6.imiss -k king.kin0 -c  MD2123_samples_change_230719.txt

# Create samples need removing

#less MD_ids_remove.txt|awk '{printf "%s\t%s\n", $1, $1}'|awk '$1 != "XXXX-I_p5" {print}' > bad_samples.txt

# keep cell line samples and corresponding MD samples so save those samples from MD_ids_remove.txt 
less MD_ids_remove_update.txt|awk '{printf "%s\t%s\n", $1, $1}' > bad_samples.txt

# Remove XXXX-I_p2 since MEGA had better genotype 
printf "XXXX-I_p2\tXXXX-I_p2\n" >> bad_samples.txt

plink \
        --bfile ${plink_prefix}_mind6 \
        --remove bad_samples.txt \
        --make-bed \
        --out ${plink_prefix}_rm_badID7

plink \
        --bfile ${plink_prefix}_rm_badID7 \
        --update-ids MD_ids_relabel.txt \
        --make-bed \
        --out ${plink_prefix}_crtID8

less ${plink_prefix}_crtID8.fam|awk '{print $1}' >MD_crtID8.list

# Samples need updating sex
Rscript 4a_MD_sr_relabel.R
awk '{print $1, $1, $2}' MD_ids_sr_relabel.txt >MD_sex_relabel.txt

# NOW samples are clean, let's update sex and check sex again, actually unnecessary to update sex since we relabelled and sex should be right.
plink \
	--bfile ${plink_prefix}_crtID8 \
	--update-sex MD_sex_relabel.txt \
	--check-sex \
	--make-bed \
	--out ${plink_prefix}_clean9

# Seperate samples by race and by mom-babies
Rscript 4b_separate_race_kin.R

# Filter by HWE by race and mom/babies
# Filter by HWE by race and mom/babies
plink \
        --bfile ${plink_prefix}_clean9 \
        --recode \
        --keep black_baby.txt \
        --hardy midp \
        --out qc_hwe_blackB

plink \
        --bfile ${plink_prefix}_clean9 \
        --recode \
        --keep black_mom.txt \
        --hardy midp \
        --out qc_hwe_blackM

plink \
        --bfile ${plink_prefix}_clean9 \
        --recode \
        --keep white_baby.txt \
        --hardy midp \
        --out qc_hwe_whiteB

plink \
        --bfile ${plink_prefix}_clean9 \
        --recode \
        --keep white_mom.txt \
        --hardy midp \
        --out qc_hwe_whiteM

# get SNPs violating HWE according to meta-p < 5*10-8 in each reported race and mom/baby

awk '{if ($9<1e-3) {print $2}}' qc_hwe_blackB.hwe > bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_whiteB.hwe >> bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_blackM.hwe >> bad_HWE_snps.list
awk '{if ($9<1e-3) {print $2}}' qc_hwe_whiteM.hwe >> bad_HWE_snps.list

# Find the overlapping SNPs appears in both mom and baby
#comm -12 <(sort bad_HWE_baby_snps.list) <(sort bad_HWE_mom_snps.list)>bad_HWE_snps.list

# Finally Carole decided to remove bad snps in either moms or babies
# filter out bad HWE SNPs
plink \
        --bfile ${plink_prefix}_clean9 \
        --exclude bad_HWE_snps.list \
        --make-bed  \
        --out ${plink_prefix}_hwe10

plink  --bfile ${plink_prefix}_hwe10 --chr 1-22 --het --out plink

Rscript --no-save 3c_heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# remove samples failed het
plink \
        --bfile ${plink_prefix}_hwe10 \
        --remove het_fail_ind.txt \
        --make-bed \
        --out ${plink_prefix}_het11


# Prune SNPs

#plink \
 #       --bfile ${plink_prefix}_het11 \
  #      --indep-pairwise 50 5 0.2 \
 #       --out indep
#
##        --split-x b37 no-fail \

#plink \
 #       --bfile ${plink_prefix}_het11 \
  #      --set-hh-missing \
   #     --make-bed \
    #    --out ${plink_prefix}_hhmissing

#plink \
 #       --bfile ${plink_prefix}_hhmissing \
  #      --extract indep.prune.in \
   #     --het \
    #    --make-bed \
     #   --out ${plink_prefix}_ld


# Estimate Ancestry PCs

# First drop ambiguous alleles from our cleaned data set and # Create a file with a list of our variants (rs#).

#plink \
 #       --bfile ${plink_prefix}_ld \
  #      --pca \
   #     --out ${plink_prefix}_pca
 
