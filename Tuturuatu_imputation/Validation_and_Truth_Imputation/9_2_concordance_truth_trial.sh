#!/bin/bash -e

# 23 April 2023
# Olivia Janes
# Imputation stats: comparing all of the trials (depth and effective popl size) within the truth imputation to the
#   chosen truth imputation (5x and 100ne) to ensure that they are all mostly the same. This ensures the chosen truth
#   imputation is the most accurate representation of the truth, regardless of depth or Ne.

#Environment: N/A  

## Edit to be run specific
sppdir=~/data/tuturuatu_all_vcf/impute/truth/
impstats=~/data/tuturuatu_all_vcf/impute/
    #Location of imputation stats folder
tlr_regions=~/data/tuturuatu_all_vcf/bcf/tlr_regions.bed
    #Define location of TLR regions bed file
run=truth

# Defining directories
impdir=${sppdir}impute/
mkdir -p ${impstats}stats/${run}_stats/concordance/
statsdir=${impstats}stats/${run}_stats/

#Ensuring gatk will run
    export PATH="~/data/programs/gatk-4.4.0.0/:$PATH"


# Concordance
    for dp in {0,4,5}
    do
        for test_ne in {default,50,100,500}
        do
            file=${impdir}beagle_imputations/filtered/Tuturuatu_VariantCalls_${dp}x_${test_ne}ne_truth_final.vcf.gz
            truth=${impdir}beagle_imputations/filtered/Tuturuatu_VariantCalls_5x_100ne_truth_final.vcf.gz
            base=$(basename ${file} _truth_final.vcf.gz)
            echo ""; echo "Calculating stats for ${base}_truth_final.vcf.gz"

            # Index according to GATK's requirements (needs a .tbi index)
            gatk IndexFeatureFile -I ${file}
            wait

            # Concordance
            gatk Concordance -eval ${file} --truth ${truth} --summary ${statsdir}concordance/${run}_${dp}x_${test_ne}ne_concordance_summary.tsv
        done
    done

# Collating concordance results
    echo ""; echo "Creating summary file for all concordance results."
    # create a new CSV file
    echo -e "Filename\ttype\tTrue Postives\tFalse Positives\tFalse Negatives\tRecall\tPrecision" > ${statsdir}concordance/truth_concordance_against_5x_100ne.csv

    # loop over all text files in the current directory
    for file in ${statsdir}concordance/*.tsv
    do
        # get the filename without the extension
        filename=$(basename ${file} _concordance_summary.tsv)
        
        # get the second row of the summary file
        row=$(sed '2!d' "${file}")
        
        # write the filename and row to the CSV file
        echo -e "${filename}\t${row}" >> ${statsdir}concordance/truth_concordance_against_5x_100ne.csv
    done

echo ""; echo "Script is complete. Find concordance results at: "
echo "${statsdir}concordance/"
echo "Find summary of concordance results at:"
echo "${statsdir}concordance/truth_concordance_against_5x_100ne.csv"
