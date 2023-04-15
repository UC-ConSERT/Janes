#!/bin/bash -e
set -e

# 13 April 2023
# Olivia Janes
# Merge the vcfs for the imputed individuals back with the high coverage individuals.

#Environment: samtools

sppdir=~/data/tuturuatu_all_vcf/
## Edit to be run specific


# Setting variables.
impdir=${sppdir}impute/
mkdir -p ${impdir}beagle_imputations/merged
impoutdir=${impdir}beagle_imputations/impute_trials/
mergedir=${impdir}beagle_imputations/merged/

for vcf in ${impoutdir}*vcf.gz
do
