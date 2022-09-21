#!/bin/bash -e 

#12 Sep 2022
#Molly Magid, adapted by Olivia Janes
#Calculating statistics to compare between filtering methods.

sppdir=~/data/tuturuatu/
statsdir=${sppdir}bcf/stats/

cd ${statsdir}

echo ", Mean Site Depth,SD,Mean Indv Depth,SD,Mean Site Missingness,SD,Mean Indv Missingness,SD,Mean Heterozygosity, Mean Heterozygosity SD">>mean_SD_filter_stats_tuturuatu_22.csv 
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
indv_miss_SD=$(awk '{x+=$5;y+=$5^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.imiss) 
 
#calculating site heterozygosity 
het_mean=$(awk '{sum +=$5} END {print sum/NR}' ${base}.het) 
het_SD=$(awk '{x+=$5;y+=$5^2}END{print sqrt(y/NR-(x/NR)^2)}' ${base}.het) 

 
#printing all stats for the file 
echo "${base}, ${site_depth_mean}, ${site_depth_SD}, ${indv_depth_mean}, \
    ${indv_depth_SD}, ${site_miss_mean}, ${site_miss_SD}, ${indv_miss_mean}, \
    ${indv_miss_SD}, ${het_mean}, ${het_SD}" >> mean_SD_filter_stats_tuturuatu_22.csv #mea$ 
done 

echo "Calculating stats and outputting into csv is done."