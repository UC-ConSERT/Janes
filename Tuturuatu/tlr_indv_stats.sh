#!/bin/bash -e

#9 Nov 2022
#Olivia Janes
#Creating individual-level stats for TLR regions within duplicate removed (nodup.bam) bam files.
#This will be used to analyse the quality of reads within the TLR regions in tuturuatu individuals, to assess if poor overall
#   genome quality correlates to poor quality within the TLR regions.
#The output file can be run through multiqc for a visual report.

## Run within the samtools environment ##

sppdir=~/data/tuturuatu/

nodupbamdir=${sppdir}nodup_bam/
        #directory that holds the merged bam files that have been sorted, fixed and had duplicates removed.
tlr_regions=${sppdir}bcf/tlr_regions.bed

mkdir -p ${nodupbamdir}tlr_indv_stats/

for file in ${nodupbamdir}*nodup.bam
do
    base=$(basename ${file} _nodup.bam)
    echo "Calculating TLR region stats and depth for ${base}..."
    samtools stats ${file} -@ 64 -t ${tlr_regions} > ${nodupbamdir}tlr_indv_stats/${base}_tlr_indv_stats.stat
    samtools depth ${file} -b ${tlr_regions} > ${nodupbamdir}tlr_indv_stats/${base}_tlr_indv_depth.stat
done

echo "Samtools stats has finished. Time to run multiqc!!!!!"


#Calculating statistics (mean depth within TLR regions) in depth.stat files.

cd ${nodupbamdir}tlr_indv_stats/

echo "Individual,Mean TLR1A Depth,SD,Mean TLR1B Depth,SD,Mean TLR2A Depth,SD,Mean TLR2B Depth,SD,Mean TLR3 Depth,SD,Mean TLR4 Depth,SD,Mean TLR5 Depth,SD,Mean TLR7 Depth,SD,Mean TLR21 Depth,SD,Mean TLR15 Depth,SD,Mean Indv Depth at TLRs,SD">>tlr_indv_depth_tuturuatu_22.csv 
for file in ${nodupbamdir}tlr_indv_stats/*indv_depth.stat
do 
    base=$(basename ${file} _tlr_indv_depth.stat)
    echo "calculating stats for ${base}" 

    #Format to follow:
    #tlr_depth_mean = sum/n
    #tlr_depth_sd = sqrt((sum^2)/n - (sum/n)^2) ##########Is this actually SD?? what is it??

    ##SD## but not actually SD
        #site_depth_SD=$(awk '{x+=$3;y+=$3^2}END{print sqrt(y/NR-(x/NR)^2)}' ${file}) 
        #y = $(awk '/jcf7180002669510/ {y+=$3^2} END {print y}' ${file})
        #x = $(awk '/jcf7180002669510/ {x+=$3} END {print x}' ${file})
        #NR = $(grep -c jcf7180002669510 ${file})
    #e.g.
        #lionsdz=$(($(awk '/lion/ {y+=$3^2} END {print y}' test1.txt)/$(grep -c lion test1.txt)-($(awk '/lion/ {x+=$3} END {print x}' test1.txt)/$(grep -c lion test1.txt))^2))
        #lionsd=$(echo "$lionsdz" | awk '{print sqrt($1)}')

    #calculating TLR1A depth 
    tlr1a_depth_mean=$(($(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})/$(grep -c jcf7180002669510 ${file})))
    #tlr1a_depth_SDz=$(($(awk '/jcf7180002669510/ {y+=$3^2} END {print y}' ${file})/$(grep -c jcf7180002669510 ${file})-($(awk '/jcf7180002669510/ {x+=$3} END {print x}' ${file})/$(grep -c jcf7180002669510 ${file}))^2))
    #tlr1a_depth_SD=$(echo "$tlr1a_depth_SDz" | awk '{print sqrt($1)}')
        ##### Issue:  not actually SD, so don't use until fixed. Leaving here so basic operation can be followed

    #calculating TLR1B depth 
    tlr1b_depth_mean=$(($(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})/$(grep -c jcf7180002669510 ${file})))

    #calculating TLR2A depth 
    tlr2a_depth_mean=$(($(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})/$(grep -c jcf7180002669510 ${file})))

    #calculating TLR1A depth 
    tlr1a_depth_mean=$(($(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})/$(grep -c jcf7180002669510 ${file})))

    #calculating TLR1A depth 
    tlr1a_depth_mean=$(($(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})/$(grep -c jcf7180002669510 ${file})))

    
    #printing all stats for the file 
    echo "${base}, ${site_depth_mean}, ${site_depth_SD}, ${indv_depth_mean}, \
        ${indv_depth_SD}, ${site_miss_mean}, ${site_miss_SD}, ${indv_miss_mean}, \
        ${indv_miss_SD}, ${het_mean}, ${het_SD}" >> mean_SD_filter_stats_tuturuatu_22.csv  
done 

echo "Calculating stats and outputting into csv is done."

#Script for running multiqc. Should be run within the multiqc environment, just copy and paste into the terminal.
#multiqc ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --outdir ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --filename multiqc_report_tlr_indv.html

#In local computer terminal, from a chosen location, use the following to extract the html report from the VM:
#rsync -rav rccuser:~/data/tuturuatu/nodup_bam/tlr_indv_stats/multiqc_report_tlr_indv.html ./