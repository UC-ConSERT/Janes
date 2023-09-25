#!/bin/bash -e 
set -e

#19 March 2023

# Olivia Janes
# A script to run all the stats for the filtered vcf from the validation run for imputation
# From: tuturuatu_imputation

# Stats scripts taken from:
#   5_0_filtering_tuturuatu_all.sh
#   5_2_snp_counts_tuturuatu_all.sh
#   5_3_calculating_stats_tuturuatu_all.sh
#   5_5_extracting_tlr_stats_tuturuatu_all.sh

## Environment: samtools

# Setting up
    run=tuturuatu_all_vcf
    ##  Needs to be edited to be run specific   ##
    sppdir=~/data/${run}/impute/validation/

    filterdir=${sppdir}bcf/filter_trial/impute/
        #directory of filtered vcf.gz

    mkdir -p ${sppdir}bcf/stats
    mkdir -p ${sppdir}bcf/stats/stats_raw_files
    statsdir=${sppdir}bcf/stats/
    cp ~/data/${run}/bcf/tlr_regions.bed ${sppdir}bcf/


#Calculating indv and site stats for pre-imputed TLR contig vcfs, to ensure that the low coverage and the validation data look similar in missingness and depth. 
echo "(5_0) Running stats script beginning."

    #Calculate stats
    echo "Calculating pre-filter stats"
    for file in ${sppdir}bcf/*concat.vcf.gz
    do
        base=$(basename ${file} _concat.vcf.gz)
        echo "Calculating depth for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base}_prefilter \
            --site-depth &
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base}_prefilter \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base}_prefilter \
            --missing-site &
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base}_prefilter \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base}_prefilter \
            --het
    done



    #calculating statistics for filtered files
    for file in ${filterdir}*.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo "Calculating depth for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base} \
            --site-depth &
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base} \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base} \
            --missing-site &
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base} \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --gzvcf ${file} \
            --out ${statsdir}stats_raw_files/${base} \
            --het
    done


echo "(5_0) Running stats script completed."
echo ""
echo ""



# (5_2) TLR SNP Counts: Compiling all SNP counts from filtering, to compare between filtering methods.
echo "(5_2) TLR SNP Counts beginning. Please fasten your seatbelts."

    bcfdir=${sppdir}bcf/
    snptxt=${statsdir}SNP_counts_validation_preimpute.txt
    tlrsnptxt=${statsdir}TLR_SNP_counts_validation_preimpute.txt
        #Defining the TLR SNP count txt output.
    #Define location of tlr_regions.bed file in script. Should be in bcf/

    ## Total SNP counts ##
    echo "Counting TOTAL SNPs"
    echo "TOTAL SNPs,Number" >> ${snptxt}
    # Prefiltered SNP counts
    #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script.
        echo "Prefiltered SNP counts"
        cd ${bcfdir}

        for file in ${bcfdir}*VariantCalls_concat.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${snptxt}
        done

    # Filtered SNP counts
        echo "Filtered SNP counts"
        cd ${filterdir}

        for file in ${filterdir}*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${snptxt}
        done


## TLR SNP counts ##
    echo "Counting TLR SNPs"
    echo "TLR SNPs,Number" >> ${tlrsnptxt}

    #Prefilter SNP counts
        echo "Prefilter TLR SNP counts"
        cd ${bcfdir}
 
        for file in ${bcfdir}*VariantCalls_concat.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${tlrsnptxt}
        done

    # Filtered SNP counts
        echo "Filtered SNP counts"
        cd ${filterdir}

        for file in ${filterdir}*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${tlrsnptxt}
        done



echo " (5_2) SNP counting is complete. Yay! Find output file at ${snptxt}"
echo ""
echo ""



#Calculating statistics to compare between filtering methods.
echo " (5_3) Calculating stats beginning. Please do not wear 3D glasses."

    statscsv=${statsdir}mean_SD_filter_stats_${run}_validation_23.csv
    ##  Needs to be edited to be run specific   ##

    cd ${statsdir}stats_raw_files/

    echo ", Mean Site Depth,SD,Mean Indv Depth,SD,Mean Site Missingness,SD,Mean Indv Missingness,SD,Mean Heterozygosity, Mean Heterozygosity SD">>${statscsv} 
    for file in ${statsdir}stats_raw_files/*.ldepth
        do 
        base=$(basename ${file} .ldepth)
        echo "calculating stats for ${base}" 
        
        #calculating site depth 
        site_depth_mean=$(awk '{sum +=$3} END {print sum/NR}' ${file}) 
        site_depth_SD=$(awk '{x+=$3;y+=$3^2}END{print sqrt(y/NR-(x/NR)^2)}' ${file}) 
        
        #calculating indv depth 
        indv_depth_mean=$(awk '{sum +=$3} END {print sum/NR}' ${base}.idepth) 
        indv_depth_SD=$(awk '{x+=$3;y+=$3^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.idepth) 
        
        #calculating site missingness 
        site_miss_mean=$(awk '{sum +=$6} END {print sum/NR}' ${base}.lmiss) 
        site_miss_SD=$(awk '{x+=$6;y+=$6^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.lmiss) 
        
        #calculating indv missingness 
        indv_miss_mean=$(awk '{sum +=$5} END {print sum/NR}' ${base}.imiss) 
        indv_miss_SD=$(awk '{x+=$5;y+=$5^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.imiss) 
        
        #calculating site heterozygosity 
        het_mean=$(awk '{sum +=$5} END {print sum/NR}' ${base}.het) 
        het_SD=$(awk '{x+=$5;y+=$5^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.het) 

        
        #printing all stats for the file 
        echo "${base}, ${site_depth_mean}, ${site_depth_SD}, ${indv_depth_mean}, \
            ${indv_depth_SD}, ${site_miss_mean}, ${site_miss_SD}, ${indv_miss_mean}, \
            ${indv_miss_SD}, ${het_mean}, ${het_SD}" >> ${statscsv}  
    done 

    echo ""
    echo "Calculating stats and outputting into csv is done."
    echo "Output file can be found at ${statsdir}${statscsv}"

echo " (5_3) Calculating stats is complete. Woopee! Find output file at ${statscsv}"
echo ""
echo ""



#Extracting the stats for TLR sites/regions from all sites stats into new files to detect any TLR sites that are outliers (poor depth/missingness/quality, etc)
echo " (5_5) Extracting TLR stats is beginning. Refrain from patting the TL-Rex for fear she loses her concentration and messes up the stats."


    mkdir -p ${statsdir}tlr_stats/

    #For prefilter stats. Will create empty files for indv stats, just ignore as we are looking for site (TLR SNP) stats
    for file in ${statsdir}stats_raw_files/*_prefilter*
    do
        echo "Finding TLR stats for ${file}"
        base=$(basename ${file})
        awk 'NR==1 {print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR1A" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002669510" && $2 > 0 && $2 < 2191) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR1B" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002696225" && $2 > 79300 && $2 < 80387) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR2A" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002687310" && $2 > 7857 && $2 < 10285) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR2B" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002686685" && $2 > 232999 && $2 < 235155) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 279461 && $2 < 279884) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 284293 && $2 < 286123) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 286306 && $2 < 286420) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR4" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002688524" && $2 > 101290 && $2 < 103549) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR5" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002688267" && $2 > 101472 && $2 < 103857) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR7" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002696332" && $2 > 371554 && $2 < 374695) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR21" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002694481" && $2 > 2000 && $2 < 3649) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR15_part" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693589" && $2 > 0 && $2 < 2103) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
    done

    ### FILTERS ###
    for file in ${statsdir}stats_raw_files/*0.6SP*
    do
        echo "Finding TLR stats for ${file}"
        base=$(basename ${file})
        awk 'NR==1 {print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR1A" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002669510" && $2 > 0 && $2 < 2191) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt

        echo "TLR1B" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002696225" && $2 > 79300 && $2 < 80387) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR2A" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002687310" && $2 > 7857 && $2 < 10285) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR2B" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002686685" && $2 > 232999 && $2 < 235155) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 279461 && $2 < 279884) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 284293 && $2 < 286123) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR3" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693511" && $2 > 286306 && $2 < 286420) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR4" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002688524" && $2 > 101290 && $2 < 103549) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR5" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002688267" && $2 > 101472 && $2 < 103857) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR7" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002696332" && $2 > 371554 && $2 < 374695) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR21" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002694481" && $2 > 2000 && $2 < 3649) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        
        echo "TLR15_part" >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
        awk '{if ($1 == "jcf7180002693589" && $2 > 0 && $2 < 2103) print $0}' ${file} >> ${statsdir}tlr_stats/${base}_tlr_regions.txt
    done

    echo "Removing (empty) individual files"
    rm ${statsdir}tlr_stats/*imiss_tlr_regions.txt
    rm ${statsdir}tlr_stats/*idepth_tlr_regions.txt

echo " (5_5) Extracting TLR stats script has finished. Now you may offer the TL-Rex some home baking or perhaps a refreshing drink."
echo ""
echo ""

echo "Entire script completed. Wow, I need a rest now!"