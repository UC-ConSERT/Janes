#!/bin/bash -e 

# 17 March 2023
# Olivia Janes

## SCRIPT DOESN'T WORK ##

# This script can be used in the home computer terminal to extract all of the outputs of the 
#   pre- ("dup") and post- ("nodup") tidyup stats runs.

homedir=~/Documents/Tuturuatu_resources/tuturuatu_all/
sppdir=rccuser:/home/rccuser/data/tuturuatu_all/    #define the spp directory path into the VM
    ## Must be edited ##

cd ${homedir}

mkdir -p nodup_bam_stats dup_bam_stats
mkdir -p dup_bam_stats/graphmap nodup_bam_stats/graphmap
mkdir -p dup_bam_stats/graphmap/coverage_histograms nodup_bam_stats/graphmap/coverage_histograms

nodupdir=${homedir}nodup_bam_stats/
dupdir=${homedir}dup_bam_stats/


##  Graphmaps   ##

#   Dup graphmaps
    echo "Extracting dup graphmaps"
    rsync -rav ${sppdir}/dup_bam_stats/*graphmap ${dupdir}graphmap/  #copy dup graphmap folders into home computer

    for file in ${dupdir}graphmap/*.graphmap    #define file as the sample's graphmap folders
    do 
        base=$(basename $file _dup.graphmap)  #pull out the sample ID to be used to rename the image
        cp $file/images_qualimapReport/*50_histogram.png ${dupdir}coverage_histograms/${base}_genome_coverage_0to50_histogram.png
            #copy the histograms and put them in the right folder, and rename
    done


#   Nodup graphmaps
    echo "Extracting nodup graphmaps"
    rsync -rav ${sppdir}/nodup_bam_stats/*graphmap ${nodupdir}graphmap/  #copy nodup graphmap folders into home computer

    for file in ${nodupdir}graphmap/*.graphmap    #define file as the sample's graphmap folders
    do 
        base=$(basename $file _nodup.graphmap)  #pull out the sample ID to be used to rename the image
        cp $file/images_qualimapReport/*50_histogram.png ${nodupdir}coverage_histograms/${base}_genome_coverage_0to50_histogram.png
            #copy the histograms and put them in the right folder, and rename
    done

echo "Script complete."