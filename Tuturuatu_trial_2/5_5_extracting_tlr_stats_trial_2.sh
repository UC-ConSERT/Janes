#!/bin/bash -e

#21 Sep 2022
#Olivia Janes
#Extracting the stats for TLR sites/regions from all sites stats into new files to detect any TLR sites that are outliers (poor depth/missingness/quality, etc)

sppdir=~/data/tuturuatu_trial_2/
statsdir=${sppdir}bcf/stats/

mkdir -p ${statsdir}tlr_stats/

#For prefilter stats. Will create empty files for indv stats, just ignore as we are looking for site (TLR SNP) stats
for file in ${statsdir}*_prefilter*
do
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

### FILTER 1 ###
#For chosen filter stats. Edit "for" statement to be chosen filter specific.
for file in ${statsdir}/Tuturuatu_VariantCalls_4x_coverage_0.1site_missing_MinGQ10.bcf.recode.*
do
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


### FILTER 2 ###
#For chosen filter stats. Edit "for" statement to be chosen filter specific.
for file in ${statsdir}/Tuturuatu_VariantCalls_4x_coverage_0.1site_missing_MinGQ20.bcf.recode.*
do
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


### FILTER 3 ###
#For chosen filter stats. Edit "for" statement to be chosen filter specific.
for file in ${statsdir}/Tuturuatu_VariantCalls_5x_coverage_0.1site_missing_MinGQ10.bcf.recode.*
do
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

echo "Script has finished"