#!/bin/bash -e 

#This file was pre changing how sed works in the for loop. This might work. use if renaming file doesn't work
#The difference is that in this file the series of sed to remove three digit numbers DOESNT change name=. The other file goes name=, name1=, name2= in a series.

#18 Nov 2022
#Olivia Janes

##  Renaming files BEFORE alignment ##

sppdir=~/data/tuturuatu_trial_2/

datadir=${sppdir}fq_gz_files/
         #directory with trimmed fastq data


#rename files to remove unnecessary text
        ##### Must be edited to be run specific #####

# I need to merge Apr and Aug wild samples
##  Apr = 20Nov19000005_A09-L1_S54_L003_R1_001.fastq.gz to  20Nov19000005_H11-L1_S125_L004_R2_001.fastq.gz
##  Aug = 20Nov19000005_A09-L1_S7_L001_R1_001.fastq.gz  to  20Nov19000005_H11-L1_S58_L001_R2_001.fastq.gz.md5


#Renaming 2019_IKMB (I164xx...) samples
#name begins as I16xx-L1_Sxxx_L003_val_[1 or 2].fq.gz
for sample in ${datadir}I*_L003_val_1.fq.gz
do 
        echo $sample
        base=$(basename $sample _val_1.fq.gz)
    #base=I164xx-L1_Sxxx_L003
        echo $base
        name=$(echo $base | sed 's/-L1_S[0-9][0-9][0-9]_L003/-L1_S00_L003/g')
        name=$(echo $base | sed 's/-L1_S[0-9][0-9]_L003//g')
    #name=I16xx
        echo $name
    #replaces ${base} with ${name}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name}/g" ${datadir}/${base}*
    #should now be I16xx_val_x.fastq.gz
done

#Renaming 2021_IKMB (Apr) samples
#names begin as 20Nov19000005_Xxx-L1_Sxxx_L00[3 or 4]_val_[1 or 2].fq.gz
for sample in ${datadir}20*L00{3..4}_val_1.fq.gz
do 
        echo $sample
        base=$(basename $sample _val_1.fq.gz)
    #base=20Nov19000005_Xxx-L1_Sxxx_L00x
        echo $base
        name1=$(echo $base | sed 's/20Nov19000005_//g')
        echo $name1
    #name1=Xxx-L1_Sxxx_L00x
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9][0-9]_L00[3-4]/-L1_S00_L003/g')
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9]_L00[3-4]/_apr/g')
    #name2=Xxx_apr
        echo $name2
    #replaces ${base} with ${name2}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name2}/g" ${datadir}/${base}*
    #should now be Xxx_apr_val_x.fastq.gz
done

#Renaming 2021_IKMB (Aug) samples
for sample in ${datadir}20*L001_val_1.fq.gz
do 
        echo $sample
        base=$(basename $sample _val_1.fq.gz)
    #base=20Nov19000005_Xxx-L1_Sxxx_L001
        echo $base
        name1=$(echo $base | sed 's/20Nov19000005_//g')
        echo $name1
    #name1=Xxx-L1_Sxxx_L001
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9][0-9]_L001/-L1_S00_L001/g')
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9]_L001/-L1_S0_L001/g')
        name2=$(echo $name1 | sed 's/-L1_S[0-9]_L001/_aug/g')
    #name2=Xxx_aug
        echo $name2
    #replaces ${base} with ${name2}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name2}/g" ${datadir}/${base}*
    #should now be Xxx_aug_val_x.fastq.gz
done

echo "Now changing val_1/2 to R1/2"

for sample in ${datadir}*.fq.gz
do
    base=$(basename $sample .fq.gz)
    name=$(echo $base | sed 's/val_/R/g')
    echo $name
    rename "s/${base}/${name}/g" ${datadir}/${base}*
done

echo "Script is complete."