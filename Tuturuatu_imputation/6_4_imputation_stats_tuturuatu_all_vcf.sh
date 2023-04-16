#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific
#r

# Defining directories
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
statsdir=${impdir}stats/
mkdir-p ${impdir}vcf_finals/vcf_merged/
mergedir=${impdir}vcf_finals/vcf_merged/


#Merge the Pre-imputation, TLR contig-seperated, study vcf back into one file for all TLR contigs
    for dp in {0,4,5}
    do
        for subset_1 in ${impdir}vcf_finals/*${dp}x_1_study.vcf.gz
        do
            base=$(basename ${subset_1} _1_study.vcf.gz)
            subset_2=${impdir}vcf_finals/${base}_2_study.vcf.gz
            subset_3=${impdir}vcf_finals/${base}_3_study.vcf.gz

            #Create a file list of files to merge
            echo ""; echo "Merging subsetted tlr contigs in:"
            echo "${subset_1}"; echo "${subset_1}" > ${mergedir}merge_list.txt
            echo "${subset_2}"; echo "${subset_2}" >> ${mergedir}merge_list.txt
            echo "${subset_3}"; echo "${subset_3}" >> ${mergedir}merge_list.txt

            #Merge and index
            bcftools concat -O z --threads 16 -f ${mergedir}merge_list.txt -o ${mergedir}${base}_study_merged.vcf.gz
            echo "Indexing ${base}_study_merged.vcf.gz"
            bcftools index -f --threads 16 ${mergedir}${base}_study_merged.vcf.gz
        done
    done


# Stats
    # To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/tuturuatu_beagle_imp.vcf.gz

    for dp in {0,4,5}
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