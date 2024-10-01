#!/bin/bash -l

#SBATCH --job-name=apc1kg
#SBATCH --nodes=1
#SBATCH --partition=tier3q
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=32gb
#SBATCH --error=/gpfs/.../log/4ancestryPC_1KG.err
#SBATCH --output=/gpfs/.../log/4ancestryPC_1KG.out

module load gcc
module load plink/1.9

# Set up
qcdir='/gpfs/.../MD2123/qcdir'
refdir='/gpfs/.../reference/reference_samples/1Kgenome/CGRCh38'
name='MD2123'
refname='all_hg38'

mkdir $qcdir/plink_log

cd ${qcdir}

# Convert variants ID to chrom:pos:A:G since coordinate ID shown in reference bim file
less $qcdir/${name}_het11.bim|awk '{print $1":"$4":"$6":"$5}' >$nameID.new
less $qcdir/${name}_het11.bim|awk '{print $2}' >$nameID.old
paste $nameID.old $nameID.new -d "\t" > $name.id

plink --bfile $qcdir/${name}_het11 \
      --update-name $name.id \
      --make-bed \
      --out $name

# Filter reference and study data for non A-T or G-C SNPs

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink --bfile  $refdir/$refname \
      --exclude $qcdir/$refname.ac_gt_snps \
      --allow-extra-chr \
      --make-bed \
      --out $qcdir/$refname.no_ac_gt_snps

mv  $qcdir/$refname.no_ac_gt_snps.log $qcdir/plink_log/$refname.no_ac_gt_snps.log

plink --bfile  $qcdir/$name \
      --exclude $qcdir/$name.ac_gt_snps \
      --make-bed \
      --out $qcdir/$name.no_ac_gt_snps

mv  $qcdir/$name.no_ac_gt_snps.log $qcdir/plink_log/$name.no_ac_gt_snps.log

# Prune study data

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --indep-pairwise 50 5 0.1 \
      --out $qcdir/$name

#mv  $qcdir/$name.prune.log $qcdir/plink_log/$name.prune.log

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --extract $qcdir/$name.prune.in \
      --make-bed \
      --out $qcdir/$name.pruned

mv  $qcdir/$name.pruned.log $qcdir/plink_log/$name.pruned.log

# Filter reference data for the same SNP set as in study
plink --bfile  $refdir/$refname \
      --extract $qcdir/$name.prune.in \
      --allow-extra-chr \
      --make-bed \
      --out $qcdir/$refname.pruned

mv  $qcdir/$refname.pruned.log $qcdir/plink_log/$refname.pruned.log

# Check and correct chromosome mismatch

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

plink --bfile $qcdir/$refname.pruned \
      --update-chr $qcdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$refname.updateChr

mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log

# Position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/${refname}.toUpdatePos


# Possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/$refname.toFlip

# Upate positions and flip alleles
plink --bfile $qcdir/$refname.updateChr \
      --update-map $qcdir/$refname.toUpdatePos 1 2 \
      --flip $qcdir/$refname.toFlip \
      --make-bed \
      --out $qcdir/$refname.flipped
mv $qcdir/$refname.flipped.log $qcdir/plink_log/$refname.flipped.log

# Remove mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
    $qcdir/$refname.mismatch

plink --bfile $qcdir/$refname.flipped \
      --exclude $qcdir/$refname.mismatch \
      --make-bed \
      --out $qcdir/$refname.clean
mv $qcdir/$refname.clean.log $qcdir/plink_log/$refname.clean.log

# Merge study genotypes and reference data
plink --bfile $qcdir/$name.pruned  \
      --bmerge $qcdir/$refname.clean.bed $qcdir/$refname.clean.bim \
         $qcdir/$refname.clean.fam  \
      --make-bed \
      --out $qcdir/$name.merge.$refname
mv $qcdir/$name.merge.$refname.log $qcdir/plink_log

# PCA on the merged data
plink --bfile $qcdir/$name.merge.$refname \
      --pca \
      --out $qcdir/$name.$refname
mv $qcdir/$name.$refname.log $qcdir/plink_log
