#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Imputing with QUILT
# Imputing tuturuatu_all bcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call bcf that has been filtered for 
#   4x coverage, 0.2 site missingness, 10 minGQ and no LD.

#Environment: stitch

sppdir=~/data/tuturuatu_all/
## Edit to be run specific

# Making a directory to hold the imputation work.
impdir=${sppdir}impute/
subsetdir=${impdir}bcf_subsets/

# Activating STITCH (required for QUILT)
R -e 'library("STITCH")'
## Is this necessary?