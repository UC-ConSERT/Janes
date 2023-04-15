#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
statsdir=${impdir}stats


# Stats
    # To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/tuturuatu_beagle_imp.vcf.gz
HAVE NOT EDITED 
    for file in ${finaldir}*_study.vcf.gz
    do
        for test_ne in {50,100,500}
        do
            base=$(basename ${file} _study.vcf.gz)
            echo ""; echo "Calculating stats for ${base}_${test_ne}ne_beagle_imp.vcf.gz"

            # Concordance
            vcf-compare ${file} ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_concordance.txt

            # Allelic/Dosage r^2
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' \
                ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_r2.txt

        done
    done

#Investigate SNPs
    #Extract information on the haplotypes at each TLR SNP for each individual
        bcftools query -R ~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' \
            ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_5x_0.6SP.vcf.gz >> ~/data/tuturuatu_all_vcf/impute/stats/tlr_haps_preimpute_5x.txt
    #Extract the headers to add to the above
        bcftools view -h ~/data/tuturuatu_all_vcf/bcf/filter_trial/impute/Tuturuatu_VariantCalls_5x_0.6SP.vcf.gz \
            | tail -n 1 >> ~/data/tuturuatu_all_vcf/impute/stats/tlr_haps_preimpute_5x_header.txt
    #Download these and extract into a spreadsheet to analyse

        bcftools query -R ~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' \
            ~/data/tuturuatu_all_vcf/impute/beagle_imputations/Tuturuatu_VariantCalls_5x_3_100ne_beagle_imp.vcf.gz >> ~/data/tuturuatu_all_vcf/impute/stats/tlr_haps_impute_5x_3.txt
    #Extract the headers to add to the above
        bcftools view -h ~/data/tuturuatu_all_vcf/impute/beagle_imputations/Tuturuatu_VariantCalls_5x_1_100ne_beagle_imp.vcf.gz \
            | tail -n 1 >> ~/data/tuturuatu_all_vcf/impute/stats/tlr_haps_impute_5x_header.txt

echo ""
echo "Imputation stats script is complete."