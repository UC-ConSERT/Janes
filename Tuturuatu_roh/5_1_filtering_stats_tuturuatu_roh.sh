#!/bin/bash -e 
set -e

#04 May 2023

# Olivia Janes
# A script to run all the stats for the filtered vcf for ROH

# Stats scripts taken from:
#   5_0_filtering_tuturuatu_all.sh
#   5_2_snp_counts_tuturuatu_all.sh
#   5_3_calculating_stats_tuturuatu_all.sh

run=tuturuatu_roh
##  Needs to be edited to be run specific   ##
sppdir=~/data/${run}/

filterdir=${sppdir}bcf/filter_trial/
    #directory of filtered vcf.gz

mkdir -p ${sppdir}bcf/stats/filter_stats/
mkdir -p ${sppdir}bcf/stats/filter_stats/stats_raw_files
statsdir=${sppdir}bcf/stats/filter_stats/


#Calculating indv and site stats for vcfs. 
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



# (5_2) SNP Counts: Compiling all SNP counts from filtering, to compare between filtering methods.
echo "(5_2) SNP Counts beginning. Please fasten your seatbelts."

    bcfdir=${sppdir}bcf/
    snptxt=${statsdir}SNP_counts_${run}.txt

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


echo " (5_2) SNP counting is complete. Yay! Find output file at ${snptxt}"
echo ""
echo ""



#Calculating statistics to compare between filtering methods.
echo " (5_3) Calculating stats beginning. Please do not wear 3D glasses."

    statscsv=${statsdir}mean_SD_filter_stats_${run}_23.csv
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


echo "Entire script completed. Wow, I need a rest now!"