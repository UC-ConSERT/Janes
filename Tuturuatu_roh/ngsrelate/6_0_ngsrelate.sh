#!/bin/bash -e 
set -e

#09 May 2023
#Olivia Janes
#Calculating inbreeding from ngsrelate using filtered vcf file.

#Environment: ngsrelate??

sppdir=~/data/tuturuatu_roh/
    ## Must be edited to be run specific
ngsrelate=~/data/programs/ngsRelate/ngsRelate

#Defining directories
    filterdir=${sppdir}bcf/filter_trial/
        #Directory containing the filtered vcf(s)
    mkdir -p ${sppdir}ngsrelate/
    ngsdir=${sppdir}ngsrelate/

#Calculating inbreeding

    for file in ${filterdir}*5x_coverage_0.2site_missing_0.6SP.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Calculating inbreeding for ${base}.vcf.gz"
        ${ngsrelate} -h ${file} -T GT -c 1 -F 1 -O ${ngsdir}${base}.res
    done

#Need to try w/o called GT: remove -T GT and -c 1
    for file in ${filterdir}*5x_coverage_0.2site_missing_0.6SP.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Calculating inbreeding for ${base}.vcf.gz"
        ${ngsrelate} -h ${file} -F 1 -O ${ngsdir}${base}_nocall.res
    done

echo "ngsRelate script has finished running."