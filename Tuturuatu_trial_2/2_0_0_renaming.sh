#!/bin/bash -e 

#18 Nov 2022
#Olivia Janes

##  Renaming files BEFORE alignment ##

sppdir=~/data/tuturuatu_trial_2/

datadir=${sppdir}fq_gz_files/
         #directory with trimmed fastq data


#rename files to remove unnecessary text
        ##### Must be edited to be run specific #####

###############to do:
#       will they be fastqs not fqs?


#Renaming 2019_IKMB (I164xx...) samples
for sample in ${datadir}I*_L003_R1_001.fq.gz
do 
        echo $sample
        base=$(basename $sample _R1_001.fq.gz)
    #base=I164xx-L1_Sxxx_L003
        echo $base
        name=$(echo $base | sed 's/-L1_S[0-9][0-9][0-9]_L003//g')
        name=$(echo $base | sed 's/-L1_S[0-9][0-9]_L003//g')
    #name=I16xx
        echo $name
    #replaces ${base} with ${name}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name}/g" ${datadir}/${base}*
    #should now be I16xx_Rx_001.fastq.gz
done

#Renaming 2021_IKMB (Apr) samples
for sample in ${datadir}20*L00{3..4}_R1_001.fq.gz
do 
        echo $sample
        base=$(basename $sample _R1_001.fq.gz)
    #base=20Nov19000005_Xxx-L1_Sxxx_L00x
        echo $base
        name1=$(echo $base | sed 's/20Nov19000005_//g')
        echo $name1
    #name1=Xxx-L1_Sxxx_L00x
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9][0-9]_L00[3-4]/_apr/g')
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9]_L00[3-4]/_apr/g')
    #name2=Xxx_apr
        echo $name2
    #replaces ${base} with ${name2}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name2}/g" ${datadir}/${base}*
    #should now be Xxx_apr_Rx_001.fastq.gz
done

#Renaming 2021_IKMB (Aug) samples
for sample in ${datadir}20*L001_R1_001.fq.gz
do 
        echo $sample
        base=$(basename $sample _R1_001.fq.gz)
    #base=20Nov19000005_Xxx-L1_Sxxx_L001
        echo $base
        name1=$(echo $base | sed 's/20Nov19000005_//g')
        echo $name1
    #name1=Xxx-L1_Sxxx_L001
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9][0-9]_L001/_aug/g')
        name2=$(echo $name1 | sed 's/-L1_S[0-9][0-9]_L001/_aug/g')
    #name2=Xxx_aug
        echo $name2
    #replaces ${base} with ${name2}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name2}/g" ${datadir}/${base}*
    #should now be Xxx_aug_Rx_001.fastq.gz
done

echo "Now removing _001"

for sample in ${datadir}*_001.fq.gz
do
    base=$(basename $sample .fq.gz)
    name=$(echo $base | sed 's/_001//g')
    echo $name
    rename "s/${base}/${name}/g" ${datadir}/${base}*
done

echo "Script is complete."