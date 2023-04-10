#!/bin/bash -e 

#29 March 2023

# Olivia Janes
# A script to download all of the stats files from the VM to the home computer
## Run this from the home computer directory for the species and run, for example: ~/Documents/Tuturuatu_resources/tuturuatu_all_rm_bad/

spprun=tuturuatu_all_vcf
##  Needs to be edited to be species and run specific

#Setting up
    cd ~/Documents/Tuturuatu_resources/
    mkdir ${spprun}
    cd ~/Documents/Tuturuatu_resources/${spprun}
    mkdir -p filter_stats
    mkdir -p filter_stats/stats_raw_files filter_stats/tlr_stats filter_stats/strand_bias_filter_stats
    mkdir -p filter_stats/strand_bias_filter_stats/tlr_stats filter_stats/strand_bias_filter_stats/sb_stats_raw_files

<<"COMMENTS"
#Download filter stats
    echo "Downloading filter stats"
    cd ~/Documents/Tuturuatu_resources/${spprun}/filter_stats/
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/stats/mean_SD_filter_stats_${spprun}_23.csv ./
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/stats/*Calls* ./stats_raw_files/
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/stats/*NP_counts.txt ./
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/stats/tlr_stats/* ./tlr_stats/
    echo ""
COMMENTS

#Download strand bias filter stats
    echo "Downloading strand bias filter stats"
    cd ~/Documents/Tuturuatu_resources/${spprun}/filter_stats/strand_bias_filter_stats/
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/strand_bias_filter_stats/mean_SD_filter_stats_${spprun}_strand_bias_23.csv ./
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/strand_bias_filter_stats/*Calls* ./sb_stats_raw_files/
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/strand_bias_filter_stats/*NP_counts_strand_bias.txt ./
    rsync -rav rccuser:/home/rccuser/data/${spprun}/bcf/strand_bias_filter_stats/tlr_stats/* ./tlr_stats/
    echo ""

echo "Script complete."