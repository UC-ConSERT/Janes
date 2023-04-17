#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Imputation stats

#Environment: samtools

## Edit to be run specific
sppdir=~/data/tuturuatu_all_vcf/
tlr_regions=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
    #Define location of TLR regions bed file
run=low_cov


# Defining directories
impdir=${sppdir}impute/
finaldir=${impdir}vcf_finals/
mkdir -p ${impdir}vcf_finals/vcf_merged/
mergedir=${impdir}vcf_finals/vcf_merged/

mkdir -p ${impdir}stats/
mkdir -p ${impdir}stats/beagle_imp_stats/
statsdir=${impdir}stats/beagle_imp_stats/



# To have a look at the imputation -> this prints it all out
    #zless -S ${impdir}beagle_imputations/filtered/[imputation.vcf.gz]


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

#Investigating TLR Haplotypes: Extract the TLR Haplotypes out of the pre-imputed and imputed files
    for dp in {0,4,5}
    do
        echo ""; echo "Extracting TLR haplotypes for ${dp}x files, pre and post impute"
        #Preimpute Haplotypes
            file=${mergedir}Tuturuatu_VariantCalls_${dp}x_study_merged.vcf.gz
            #Extract information on the haplotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file} > ${statsdir}tlr_haps_preimpute_${dp}x_${run}.txt
            #Extract the headers to add to the above
                bcftools view -h ${file} | tail -n 1 > ${statsdir}tlr_haps_preimpute_header_${run}.txt
            #Download these and extract into a spreadsheet to analyse

        #Imputed Haplotypes (for 100ne only)
            file2=${impdir}beagle_imputations/filtered/Tuturuatu_VariantCalls_${dp}x_100ne_filtered.vcf.gz
            #Extract information on the haplotypes at each TLR SNP for each individual
                bcftools query -R ${tlr_regions} --format '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' ${file2} > ${statsdir}tlr_haps_impute_${dp}x_${run}.txt
            #Extract the headers to add to the above
                bcftools view -h ${file2} | tail -n 1 > ${statsdir}tlr_haps_impute_header_${run}.txt
            #Download these and extract into a spreadsheet to analyse
    done

# Stats
    
NOT EDITED BELOW
    for dp in {0,4,5}
    do
        for test_ne in {50,100,500}
        do
            file=${impdir}beagle_imputations/filtered/Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_filtered.vcf.gz
            base=$(basename ${file} _filtered.vcf.gz)
            echo ""; echo "Calculating stats for ${base}_filtered.vcf.gz"

            # Concordance
            vcf-compare ${file} ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_concordance.txt

            # Allelic/Dosage r^2
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\t%IMP\n' \
                ${impdir}beagle_imputations/${base}_${test_ne}ne_beagle_imp.vcf.gz > ${impdir}stats/${base}_${test_ne}ne_beagle_imp_r2.txt

        done
    done


echo ""
echo "Imputation stats script is complete."