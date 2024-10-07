# Genotype_QC

Genotypes were lifted over as necessary

Then genotypes of different arrays were first QC’d separately (by sex or kinship)

Do ancestry PCs before or after merging

Fix reference and add annotation (DB SNP) using bcftools and merge genotypes of 
different arrays

Conduct genotype QC again,  ancestry PCs and check Mendelian error as necessary

Prepare input for imputation

Do imputation on TopMed Server
