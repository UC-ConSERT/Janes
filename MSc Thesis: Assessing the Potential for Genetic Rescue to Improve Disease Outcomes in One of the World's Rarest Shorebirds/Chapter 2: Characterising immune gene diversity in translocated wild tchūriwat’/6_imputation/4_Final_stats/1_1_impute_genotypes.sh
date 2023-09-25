#!/bin/bash -e 
set -e

# 09 May 2023

# Olivia Janes
# Investigating TLR Genotypes: Extract the TLR Genotypes
# From: tuturuatu_imputation

## Environment: samtools

# Setting up
    ##  Needs to be edited to be run specific   ##
        sppdir=~/data/tuturuatu_all_vcf/
        run=low_cov
        vcf=Tuturuatu_VariantCalls_0.9miss_low_cov_final.vcf.gz
        tlr_regions=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
            #Define location of TLR regions bed file

    #Setting directories
        statsdir=${sppdir}final_outputs/final_stats/
        workdir=${sppdir}final_outputs/final_vcfs/


#Investigating TLR Genotypes: Extract the TLR Genotypes out of the imputed at Ne=100 files
    for file in ${workdir}${vcf}
    do
        echo ""; echo "Extracting TLR Genotypes for ${dp}x files, post impute"

            #Extract information on the Genotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file} > ${statsdir}tlr_genotypes_impute_0.9miss.txt
            #Extract the headers to add to the above
                bcftools view -h ${file} | tail -n 1 > ${statsdir}tlr_genotypes_impute_header_0.9miss.txt
            #Download these and extract into a spreadsheet to analyse

    done
    echo "TLR Genotype files can be found at: ${statsdir}tlr_genotypes..."



echo ""; echo "To download all of the stats, navigate to the right directory on your desktop: ~/Documents/Tuturuatu_resources/tuturuatu_all_vcf/final_outputs/final_stats/"
echo "Enter code (edited for the right run):"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/tlr_geno* ./"

echo ""
echo "Imputation genotype script is complete."       
