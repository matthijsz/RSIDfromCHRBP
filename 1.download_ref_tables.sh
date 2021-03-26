#!/bin/bash

# Change this link depending on what GRCh/SNP build you want
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
      vcftools --gzvcf All_20180418.vcf.gz --chr $CHR --get-INFO RS --out chr$CHR
      vcftools --gzvcf All_20180418.vcf.gz --chr $CHR --get-INFO dbSNPBuildID --out buildchr$CHR
done