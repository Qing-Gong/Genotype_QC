#!/bin/bash -l

#SBATCH --job-name=Lift-over
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=32gb
#SBATCH --output=/scratch/qgong/.../log/2MP_liftOver.out
#SBATCH --error=/scratch/qgong/.../log/2MP_liftOver.err

#Set arguments
plink_prefix=MP_updSex1 
out_dir=/scratch/qgong/.../getPlink

cd ${out_dir}

module load gcc/11.3.0
module load crossmap
module load R
module load plink/1.9
module load king

#Remove SNPs with duplicate positions

cat MD2123.dupvar | sed -e '1d' | cut -f4 > tmp_dupl_snpids.txt
plink \
	--bfile $plink_prefix \
	--not-chr 0, 25-26 \
	--exclude tmp_dupl_snpids.txt \
        --make-bed \
	--out ${plink_prefix}_0dupl

#Create bed file to crossover from hg19 to hg38 
cat ${plink_prefix}_0dupl.bim | cut -f1 |sed 's/23/X/'|sed 's/24/Y/'| sed 's/^/chr/'> tmp_c1.txt
cat ${plink_prefix}_0dupl.bim | cut -f4 > tmp_c2.txt
cat ${plink_prefix}_0dupl.bim | cut -f4 > tmp_c3.txt
cat ${plink_prefix}_0dupl.bim | cut -f2 > tmp_c4.txt
paste  tmp_c1.txt \
       tmp_c2.txt \
       tmp_c3.txt \
       tmp_c4.txt \
       >  tmp_in.bed

#Do crossover
CrossMap.py bed hg19ToHg38.over.chain \
            tmp_in.bed  \
            tmp_out.bed

#Extract only those SNPs that were successfully cross-overed
cut -f4 tmp_out.bed > tmp_snp_keep.txt
plink \
	--bfile ${plink_prefix}_0dupl \
      	--extract tmp_snp_keep.txt \
      	--make-bed \
     	 --out ${plink_prefix}_liftOver

#Update bim file positions
cat /scratch/qgong/.../update_pos.R | R --vanilla

#Remove strand ambiguous SNPs
cat /scratch/qgong/.../get_strand_amb_SNPs.R | R --vanilla

plink \
	--bfile ${plink_prefix}_liftOver \
	--exclude tmp_strand_remove_snps.txt \
	--sort-vars \
	--missing \
	--make-bed \
	--out ${plink_prefix}_liftOver_no_AT_CG

#Cleanup
rm tmp_*


king -b ${plink_prefix}_liftOver_no_AT_CG.bed --kinship
