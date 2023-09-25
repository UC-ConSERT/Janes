#!/bin/bash -e
set -e

# 23 April 2023

# Olivia Janes
# Imputation stats: comparing all of the trials (depth and effective popl size) within the downsampling imputation to the
#   chosen truth imputation (5x and 100ne) to assess the effect of downsampling.
#   Will focus on 5x, 100Ne trials, but the other depths are important for their effect when coverage is lower. The best filter might change!
# From: tuturuatu_imputation

## Environment: This MUST be run outside of the conda environment.

# Setting up
    ## Edit to be run specific
    sppdir=~/data/tuturuatu_all_vcf/
    impstats=${sppdir}impute/stats/
        #Location of imputation stats folder
    tlr_regions=${sppdir}bcf/tlr_regions.bed
        #Define location of TLR regions bed file
    truthvcf=${sppdir}impute/truth/impute/beagle_imputations/filtered/Tuturuatu_VariantCalls_5x_100ne_truth_final.vcf.gz
        #Chosen truth vcf
    run=downsampling


for ds in {0.05,0.1,0.2,0.3}
do

    echo "Beginning script for ${ds}"

    
        evaldir=${sppdir}impute/downsampling_trial/${ds}_downsample/impute/beagle_imputations/filtered/
            #Directory of eval/test vcfs to be tested against the truth vcf
        mkdir -p ${impstats}${run}_stats/${ds}_downsample_stats/concordance/
        statsdir=${impstats}${run}_stats/${ds}_downsample_stats/
        summarycsv=${statsdir}concordance/${run}_${ds}_concordance_against_5x_100ne_truth.csv
            #Csv to hold the collated summary stats for the concordance analysis

    #Ensuring gatk will run
        export PATH="~/data/programs/gatk-4.4.0.0/:$PATH"


    # Concordance
        for dp in {0,4,5}
        do
            for test_ne in {default,50,100,500}
            do
                file=${evaldir}Tuturuatu_VariantCalls_${ds}ds_${dp}x_${test_ne}ne_${run}_final.vcf.gz
                base=$(basename ${file} _${run}_final.vcf.gz)
                echo ""; echo "Calculating stats for ${base}_${run}_final.vcf.gz"

                # Index according to GATK's requirements (needs a .tbi index)
                gatk IndexFeatureFile -I ${file}
                wait

                # Concordance
                gatk Concordance -eval ${file} --truth ${truthvcf} --summary ${statsdir}concordance/${run}_${ds}_${dp}x_${test_ne}ne_concordance_summary.tsv
            done
        done

    # Collating concordance results
        echo ""; echo "Creating summary file for all concordance results."
        # create a new CSV file
        echo -e "Filename\ttype\tTrue Postives\tFalse Positives\tFalse Negatives\tRecall\tPrecision" > ${summarycsv}

        # loop over all text files in the current directory
        for file in ${statsdir}concordance/*.tsv
        do
            # get the filename without the extension
            filename=$(basename ${file} _concordance_summary.tsv)
            
            # get the second row of the summary file
            row=$(sed '2!d' "${file}")
            
            # write the filename and row to the CSV file
            echo -e "${filename}\t${row}" >> ${summarycsv}
        done

done

echo ""; echo "Script is complete. Find concordance results at: "
echo "${statsdir}concordance/"
echo "Find summary of concordance results at:"
echo "${summarycsv}"
