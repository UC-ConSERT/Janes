#!/bin/bash -e 

# 16 Mar 2023
# Olivia Janes
# From: tuturuatu_all
# Creating a README.md file that describes all of the bam files in the to_merge/ folder. 

## Environment: N/A

sppdir=~/data/tuturuatu_all/

datadir=${sppdir}to_merge_bam/
         #directory with trimmed fastq data
species="Tuturuatu"

echo "Beginning the task..."

# Create a README to describe all of the files
    # Create the header of the README.md file
    echo "# File Descriptions" > ${datadir}README.md
    echo "" >> ${datadir}README.md

    # Loop through all the BAM files in the folder
    for file in ${datadir}CT{01..12}.bam
    do
    # Get the file name without the extension
    name=$(basename "$file" | cut -f 1 -d '.')

    # Get the SM tag from the header of the BAM file
    sm_tag=$(samtools view -H "$file" | grep -m 1 '^@RG' | sed -n 's/^.*SM:\([^\t]*\).*$/\1/p')

    # Add a description for the BAM file
    if [ -n "$sm_tag" ]; then
        description="BAM file for sample/individual $name. Sample ID: $sm_tag."
    else
        description="BAM file for sample/individual $name. No SM tag (Sample ID) in header."
    fi

    # Add the file name and description to the README.md file
    echo "## $name" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    echo "$description" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    done

# For m.bams

    # Loop through all the BAM files in the folder
    for file in ${datadir}*m.bam
    do
    # Get the file name without the extension
    name=$(basename "$file" | cut -f 1 -d '.')
    indv=$(basename $file m.bam)

    # Get the SM tag from the header of the BAM file
    sm_tag=$(samtools view -H "$file" | grep -m 1 '^@RG' | sed -n 's/^.*SM:\([^\t]*\).*$/\1/p')

    # Add a description for the BAM file
    if [ -n "$sm_tag" ]; then
        description="Second BAM file for sample/individual $indv. Sample ID: $sm_tag. This file is to be reheadered and then merged with $indv.bam"
    else
        description="Second BAM file for sample/individual $indv. No SM tag (Sample ID) in header. This file is to be reheadered and then merged with $indv.bam"
    fi

    # Add the file name and description to the README.md file
    echo "## $name" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    echo "$description" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    done

# For mr.bams

    # Loop through all the BAM files in the folder
    for file in ${datadir}*mr.bam
    do
    # Get the file name without the extension
    name=$(basename "$file" | cut -f 1 -d '.')
    indv=$(basename $file mr.bam)

    # Get the SM tag from the header of the BAM file
    sm_tag=$(samtools view -H "$file" | grep -m 1 '^@RG' | sed -n 's/^.*SM:\([^\t]*\).*$/\1/p')

    # Add a description for the BAM file
    if [ -n "$sm_tag" ]; then
        description="(Reheadered) second BAM file for sample/individual $indv. Sample ID: $sm_tag. This file is to be merged with $indv.bam"
    else
        description="(Reheadered) second BAM file for sample/individual $indv. No SM tag (Sample ID) in header. This file is to be merged with $indv.bam"
    fi

    # Add the file name and description to the README.md file
    echo "## $name" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    echo "$description" >> ${datadir}README.md
    echo "" >> ${datadir}README.md
    done

    echo "Complete! That definitely wasn't a waste of time!"
    echo ""
    echo "Please check ${datadir}README.md for a description of the files in ${datadir}."