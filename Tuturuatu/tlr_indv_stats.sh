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

<<"COMMENTS"

mkdir -p ${nodupbamdir}tlr_indv_stats/

for file in ${nodupbamdir}*nodup.bam
do
    base=$(basename ${file} _nodup.bam)
    echo "Calculating TLR region stats and depth for ${base}..."
    samtools stats ${file} -@ 64 -t ${tlr_regions} > ${nodupbamdir}tlr_indv_stats/${base}_tlr_indv_stats.stat
    samtools depth ${file} -b ${tlr_regions} > ${nodupbamdir}tlr_indv_stats/${base}_tlr_indv_depth.stat
done

echo "Samtools stats has finished. Time to run stats, stat!!!!!"

COMMENTS

#Calculating statistics (mean depth within TLR regions) in depth.stat files.
    #Depsite this searching the entire chromosome, these are based on stats files that only contain TLR regions, so is actually only calculating depth based on TLR regions.

cd ${nodupbamdir}tlr_indv_stats/

echo "Individual,Mean TLR1A Depth,Mean TLR1B Depth,Mean TLR2A Depth,Mean TLR2B Depth,Mean TLR3 Depth,Mean TLR4 Depth,Mean TLR5 Depth,Mean TLR7 Depth,Mean TLR21 Depth,Mean TLR15 (part) Depth,Mean Indv Depth at all TLRs">> tlr_indv_depth_stats_tuturuatu_22.csv 
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
    x_1a=$(awk '/jcf7180002669510/ {sum +=$3} END {print sum}' ${file})
    n_1a=$(grep -c jcf7180002669510 ${file})
    tlr1a_depth_mean=$(echo "$x_1a,$n_1a" | awk '{print $1/$2}')
    #tlr1a_depth_mean=$(echo "$x_1a,$n_1a" | awk 'BEGIN {print $1/$2}')
    #tlr1a_depth_SDz=$(($(awk '/jcf7180002669510/ {y+=$3^2} END {print y}' ${file})/$(grep -c jcf7180002669510 ${file})-($(awk '/jcf7180002669510/ {x+=$3} END {print x}' ${file})/$(grep -c jcf7180002669510 ${file}))^2))
    #tlr1a_depth_SD=$(echo "$tlr1a_depth_SDz" | awk '{print sqrt($1)}')
        ##### Issue:  not actually SD, so don't use until fixed. Leaving here so basic operation can be followed

    #calculating TLR1B depth 
    x_1b=$(awk '/jcf7180002696225/ {sum +=$3} END {print sum}' ${file})
    n_1b=$(grep -c jcf7180002696225 ${file})
    tlr1b_depth_mean=$(echo "$x_1b,$n_1b" | awk '{print $1/$2}')

    #calculating TLR2A depth 
    x_2a=$(awk '/jcf7180002687310/ {sum +=$3} END {print sum}' ${file})
    n_2a=$(grep -c jcf7180002687310 ${file})
    tlr2a_depth_mean=$(echo "$x_2a,$n_2a" | awk '{print $1/$2}')

    #calculating TLR2B depth 
    x_2b=$(awk '/jcf7180002686685/ {sum +=$3} END {print sum}' ${file})
    n_2b=$(grep -c jcf7180002686685 ${file})
    tlr2b_depth_mean=$(echo "$x_2b,$n_2b" | awk '{print $1/$2}')

    #calculating TLR3 depth 
    x_3=$(awk '/jcf7180002693511/ {sum +=$3} END {print sum}' ${file})
    n_3=$(grep -c jcf7180002693511 ${file})
    tlr3_depth_mean=$(echo "$x_3,$n_3" | awk '{print $1/$2}')

    #calculating TLR4 depth 
    x_4=$(awk '/jcf7180002688524/ {sum +=$3} END {print sum}' ${file})
    n_4=$(grep -c jcf7180002688524 ${file})
    tlr4_depth_mean=$(echo "$x_4,$n_4" | awk '{print $1/$2}')

    #calculating TLR5 depth 
    x_5=$(awk '/jcf7180002688267/ {sum +=$3} END {print sum}' ${file})
    n_5=$(grep -c jcf7180002688267 ${file})
    tlr5_depth_mean=$(echo "$x_5,$n_5" | awk '{print $1/$2}')

    #calculating TLR7 depth 
    x_7=$(awk '/jcf7180002696332/ {sum +=$3} END {print sum}' ${file})
    n_7=$(grep -c jcf7180002696332 ${file})
    tlr7_depth_mean=$(echo "$x_7,$n_7" | awk '{print $1/$2}')

    #calculating TLR21 depth 
    x_21=$(awk '/jcf7180002694481/ {sum +=$3} END {print sum}' ${file})
    n_21=$(grep -c jcf7180002694481 ${file})
    tlr21_depth_mean=$(echo "$x_21,$n_21" | awk '{print $1/$2}')

    #calculating TLR15 (part) depth 
    x_15=$(awk '/jcf7180002693589/ {sum +=$3} END {print sum}' ${file})
    n_15=$(grep -c jcf7180002693589 ${file})
    tlr15_depth_mean=$(echo "$x_15,$n_15" | awk '{print $1/$2}')

    #calculating TLR depth overall
    tlr_overall_depth_mean=$(awk '{sum +=$3} END {print sum/NR}' ${file}) 
    
    #printing all stats for the file 
    echo "${base}, ${tlr1a_depth_mean}, ${tlr1b_depth_mean}, ${tlr2a_depth_mean},\
        ${tlr2b_depth_mean}, ${tlr3_depth_mean}, ${tlr4_depth_mean}, ${tlr5_depth_mean}, \
        ${tlr7_depth_mean}, ${tlr21_depth_mean}, ${tlr15_depth_mean}, ${tlr_overall_depth_mean}" >> tlr_indv_depth_stats_tuturuatu_22.csv  
done 

echo "Calculating stats and outputting into csv is done."
echo "Time to run multiqc and export the stats files.
    "
echo "Running MultiQC should be done within a multiqc environment. The code is: 
    multiqc ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --outdir ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --filename multiqc_report_tlr_indv.html
    "
echo " To export the multiqc report: In local computer terminal, from a chosen location, use the following to extract the multiqc html report from the VM:
    rsync -rav rccuser:~/data/tuturuatu/nodup_bam/tlr_indv_stats/multiqc_report_tlr_indv.html ./
    "
echo "To export the depth stats: In local computer terminal, from a chosen location, use the following to extract the TLR indv depth stats csv from the VM:
    rsync -rav rccuser:~/data/tuturuatu/nodup_bam/tlr_indv_stats/tlr_indv_depth_stats_tuturuatu_22.csv ./"

## MultiQC ##
    #Script for running multiqc. Should be run within the multiqc environment, just copy and paste into the terminal.
    #multiqc ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --outdir ~/data/tuturuatu/nodup_bam/tlr_indv_stats/ --filename multiqc_report_tlr_indv.html

    #In local computer terminal, from a chosen location, use the following to extract the multiqc html report from the VM:
    #rsync -rav rccuser:~/data/tuturuatu/nodup_bam/tlr_indv_stats/multiqc_report_tlr_indv.html ./

## Depth stats ##
    #In local computer terminal, from a chosen location, use the following to extract the TLR indv depth stats csv from the VM:
    #rsync -rav rccuser:~/data/tuturuatu/nodup_bam/tlr_indv_stats/tlr_indv_depth_stats_tuturuatu_22.csv ./