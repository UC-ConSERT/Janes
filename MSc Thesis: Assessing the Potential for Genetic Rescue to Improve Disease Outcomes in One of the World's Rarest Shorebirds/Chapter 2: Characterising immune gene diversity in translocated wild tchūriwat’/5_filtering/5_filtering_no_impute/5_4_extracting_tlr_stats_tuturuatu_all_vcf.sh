#!/bin/bash -e

# 21 Sep 2022

# Olivia Janes
# Extracting the stats for TLR sites/regions from all sites stats into new files to detect any TLR sites that are outliers (poor depth/missingness/quality, etc)
# From: tuturuatu_all_vcf

## Environment: N/A

# Setting up
    sppdir=~/data/tuturuatu_all_vcf/
    statsdir=${sppdir}bcf/stats/

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
    for file in ${statsdir}stats_raw_files/*missing*
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

echo "Script has finished"