#!/bin/bash -e 
set -e

#04 May 2023
#Olivia Janes
#Calculating Runs of Homozygosity (ROH) using filtered vcf file.

#Environment: plink

sppdir=~/data/tuturuatu_roh/
    ## Must be edited to be run specific

#Defining directories
    vcf_out=${bcfdir}filter_trial/
        #Directory containing the filtered vcf(s)




echo "Script has finished running. Now whether it worked or not is another question..."