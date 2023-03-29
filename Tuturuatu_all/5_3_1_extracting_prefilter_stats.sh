#!/bin/bash -e 

#23 March 2023
#Olivia Janes
#Extracting statistics from previously made Prefilter stats files.

## This script DOESN'T work yet?


sppdir=~/data/tuturuatu_all/
statsdir=${sppdir}bcf/stats/
statscsv=${statsdir}prefilter_stats_tuturuatu_all_23.csv
##  Needs to be edited to be run specific   ##


#Setting csv headers
echo "Name,INDV,N Sites,Mean Indv Depth,N Data,N Missingness,Mean Indv Missingness">>${statscsv} 

#Getting individual names
for bam in ${sppdir}nodup_bam/*_nodup.bam
do
    name=$(basename ${bam} _nodup.bam)
    echo "${name}">>${statsdir}indv_names.txt
done

<<"COMMENTS"

for file in ${statsdir}*prefilter.idepth
do 
    base=$(basename ${file} .idepth)
    echo "extracting stats for ${file}" 

    #extracting indv depth stats
    INDV=$(awk '{x+=$1} END {print x}' ${file})
    n_sites=$(awk '{x+=$2} END {print x}' ${file})
    indv_depth_mean=$(awk '{x+=$3} END {print x}' ${file})

    #extracting indv missingness stats
    n_data=$(awk '{x+=$2} END {print $2}' ${statsdir}${base}.imiss)
    n_miss=$(awk '{x+=$4} END {print x}' ${statsdir}${base}.imiss)
    indv_miss_mean=$(awk '{x+=$5} END {print x}' ${statsdir}${base}.imiss)

    #printing all stats for the file 
    echo ",${INDV},${n_sites},${indv_depth_mean},\
        ${n_data},${n_miss},${indv_miss_mean}" >> ${statscsv}  
done 
COMMENTS

echo ""
echo "Calculating stats and outputting into csv is done."
echo "Output file can be found at ${statscsv}"