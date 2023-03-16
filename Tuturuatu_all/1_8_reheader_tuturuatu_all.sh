#!/bin/bash -e 

#16 Mar 2023
#  Reheading samples that are to be merged, to ensure that the merging is successful.
## For CAPTIVE samples from 2019 (IKMB) and 2021 (LIC) that need to be merged before variant calling with the rest. ##
## This script must be run before merging. 

#Olivia Janes

sppdir=~/data/tuturuatu_all/

datadir=${sppdir}to_merge_bam/
         #directory with trimmed fastq data
species="Tuturuatu"

## ENV must be samtools

<<"COMMENTS"

#rename files to remove unnecessary text
        ##### Must be edited to be sample specific #####
for sample in ${datadir}*_merged.bam
do 
        echo $sample
        base=$(basename $sample .bam)
        echo $base
        name=$(echo $base | sed 's/_merged//g')
        echo $name
        rename "s/${base}/${name}/g" ${datadir}/${base}* 
done

#Here I manually changed each of the I164... file names to their corresponding CT.. file name plus "m" for merging

COMMENTS

## Reheading
for file in ${datadir}*m.bam    #for files to be merged, where the header SM is the old sample ID (e.g. I164xxx)
do
    base=$(basename $file m.bam)    # get new header
    echo "Reheading $file"

    # extract old sample ID from header of the "m.bam" file. Beware, this will replace the SM tag in all of the @RG lines!
    old_id=$(samtools view -H "$file" | awk '/^@RG/ {for(i=1;i<=NF;i++) if($i~/^SM:/) {split($i,a,":"); print a[2];}}') 

    samtools view -H $file | sed "s/SM:${old_id}/SM:${base}/" \
        | samtools reheader - ${datadir}${base}m.bam > ${datadir}${base}mr.bam    # reheader the file xxxxm.bam with new header to match the other file (xxxx.bam), output as xxxxmr.bam
done

# Creating a file to double check that all SM tags are correct.
echo "Reheading complete. Creating a file of SM tags."

# Create a new text file to store the SM tags
echo -e "File Name\tSM Tag" > ${datadir}sm_tags.txt

# Loop through all the BAM files in the folder
for file in ${datadir}*.bam
do
  # Extract the SM tag from the BAM file header
  sm_tag=$(samtools view -H "$file" | grep '^@RG' | sed -n 's/^.*SM:\([^\t]*\).*$/\1/pg' | tr '\n' ',')

  # Append the file name and SM tag to the text file
  echo -e "$file\t$sm_tag" >> ${datadir}sm_tags.txt
done

echo "Script is complete. Files are ready to be merged."
echo ""
echo "Please check ${datadir}sm_tags.txt to ensure that all SM tags match the sample/indv ID."
echo "You may now delete CTxxm.bam, as CTxxmr.bam are the new files to be merged with CTxx.bam"