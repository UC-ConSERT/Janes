#!/bin/bash -e 

#19 March 2023

# Olivia Janes
# A script to run all the stats for the strand bias filtered bcfs.
# I will do this for BOTH of the methods of strand bias filtering: bcftools +setGT and bcftools query.

# Stats scripts taken from:
#   5_0_filtering_tuturuatu_all.sh
#   5_2_snp_counts_tuturuatu_all.sh
#   5_3_calculating_stats_tuturuatu_all.sh
#   5_5_extracting_tlr_stats_tuturuatu_all.sh

sppdir=~/data/tuturuatu_all/

sbiasdir=${sppdir}bcf/filter_strand_bias/
    #directory of strand bias filtered bcfs from bcftools +setGT
bquerydir=${sppdir}bcf/filter_strand_bias_bcftools_query/
    #directory of strand bias filtered bcfs from bcftools query

mkdir -p ${sppdir}bcf/strand_bias_filter_stats
statsdir=${sppdir}bcf/strand_bias_filter_stats/


echo "(5_0) Running stats script beginning."

    #calculating statistics for 'bcftools +setGT' strand bias filtered files
    for file in ${sbiasdir}*.bcf
    do
        base=$(basename ${file} .bcf)
        echo "Calculating depth for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_setgt \
            --site-depth &
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_setgt \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_setgt \
            --missing-site &
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_setgt \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_setgt \
            --het
    done

    #calculating statistics for 'bcftools query' strand bias filtered files
    for file in ${bquerydir}*.bcf
    do
        base=$(basename ${file} .bcf)
        echo "Calculating depth for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_bquery \
            --site-depth &
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_bquery \
            --depth &
        echo "Calculating missingness for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_bquery \
            --missing-site &
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_bquery \
            --missing-indv &
        echo "Calculating individual heterozygosity for ${base}..."
        vcftools --bcf ${file} \
            --out ${statsdir}${base}_bquery \
            --het
    done

echo "(5_0) Running stats script completed."
echo ""
echo ""


# (5_2) SNP Counts: Compiling all SNP counts from filtering, to compare between filtering methods.
echo "(5_2) SNP Counts beginning. Please fasten your seatbelts."

    bcfdir=${sppdir}bcf/
    snptxt=${statsdir}TLR_SNP_counts_strand_bias.txt
        #Defining the SNP count txt output.
    #Define location of tlr_regions.bed file in script. Should be in bcf/


    echo "##### Before filtering SNP count" >> ${snptxt}

    #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
    for file in ${bcfdir}*VariantCalls_concat.bcf
    do
        base=$(basename ${file} .bcf)
        echo ${base} >> ${snptxt}
        bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${snptxt}
    done


    echo "##### bcftools +setGT Filtered SNP counts" >> ${snptxt}

    for file in ${sbiasdir}/*.bcf
    do
        base=$(basename ${file} .bcf)
        bcftools index ${base}.bcf
        echo ${base} >> ${snptxt}
        bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${snptxt}
    done


    echo "##### bcftools query Filtered SNP counts" >> ${snptxt}

    for file in ${bquerydir}/*.bcf
    do
        base=$(basename ${file} .bcf)
        bcftools index ${base}.bcf
        echo ${base} >> ${snptxt}
        bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${base}.bcf | wc -l >> ${snptxt}
    done

echo " (5_2) SNP counting is complete. Yay! Find output file at ${snptxt}"
echo ""
echo ""



#Calculating statistics to compare between filtering methods.
echo " (5_3) Calculating stats beginning. Please do not wear 3D glasses."

    statscsv=${statsdir}mean_SD_filter_stats_tuturuatu_all_strand_bias_23.csv
    ##  Needs to be edited to be run specific   ##


    echo ", Mean Site Depth,SD,Mean Indv Depth,SD,Mean Site Missingness,SD,Mean Indv Missingness,SD,Mean Heterozygosity, Mean Heterozygosity SD">>${statscsv} 
    for file in ${statsdir}*.ldepth
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
    indv_miss_mean=$(awk '{sum +=$4} END {print sum/NR}' ${base}.imiss) 
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
    for file in ${sppdir}bcf/stats/*_prefilter*
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
    for file in ${statsdir}*.bcf*
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