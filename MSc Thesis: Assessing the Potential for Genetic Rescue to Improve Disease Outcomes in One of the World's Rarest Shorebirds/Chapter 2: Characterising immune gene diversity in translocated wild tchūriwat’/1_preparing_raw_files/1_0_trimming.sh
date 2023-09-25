#!/bin/bash -e

# 18 Nov 2022
# Molly Magid, small modifications made by Olivia Janes
# Script for trimming and processing 'tchuriwat/tūturuatu raw reads (Molly's 2019 Captive & 2021 Wild (IKMB))
# From: tuturuatu_trial_2

## Environment: fastqc

#directory with raw fastq data      ##### Must be edited to be run specific.
rawdata=~/data/tuturuatu_trial_2/raw_fq_gz_files/
#directory with fastq data
datadir=~/data/tuturuatu_trial_2/fq_gz_files/

#run Fastqc v. 0.11.8 on raw data files
for file in ${rawdata}*.gz
do
    fastqc $file
done

echo "FastQC complete. Now trimming..."

#trim files with TrimGalore v. 0.6.6 before alignment- change parameters based on results
for file in ${rawdata}*_R1_001.fastq.gz
do
       base=$(basename ${file} _R1_001.fastq.gz)
       trim_galore \
       --paired --nextera --cores 6 --nextseq 28 \
       --length 50 --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
       --clip_R1 20 --clip_R2 20  --2colour 20 --fastqc ${file} \
       ${rawdata}${base}_R2_001.fastq.gz -o ${datadir} --basename ${base}
       #trims paired end reads of a sample to a minimum length of 50, does a 3' clip of 5 bp and 5' clip of 20 bp,
       #performs clipps with a 2-colour compatible quality phred score of 20, and clip of 20 bases
       #also runs fastqc on trimmed files
done
wait
echo "Trimming is complete. Now compiling reports..."

#collects trimming reports into one file
for file in ${datadir}*.txt
do
cat $file >> trimming_reports.txt 
done

echo "Script is done."