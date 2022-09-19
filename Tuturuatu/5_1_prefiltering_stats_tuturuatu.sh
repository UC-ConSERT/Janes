#!/bin/bash -e 

#12 Sep 2022
#Adapted by Olivia Janes from https://speciationgenomics.github.io/filtering_vcfs/
#Pre filtering statistics calculated to guide vcf filtering and compare to post-filtering statistics.

sppdir=~/data/tuturuatu/
bcfdir=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

mkdir -p ${bcfdir}stats/

varcallbcf=${bcfdir}Tuturuatu_VariantCalls_concat.bcf
    #output concatenated file from variant calling
    ##### Must be edited to be sample specific #####
    ##### Must change to bcf or vcf as previous file was #####
outfile=${bcfdir}stats/Tuturuatu_VariantCalls_prefilter
    #defining the output file prefix

#If uncompressed vcf, should be compressed here using bgzip [variantcalls].vcf
#Index bcf
echo "Indexing variant calls file"
bcftools index ${varcallbcf}
wait

#Calculate stats
echo "Calculating pre-filter stats"

#Calculate allele frequency for each variant
vcftools --bcf ${varcallbcf} --freq2 --out ${outfile} --max-alleles 2
    ##### output may be changed to --gzvcf so check ######
    ##### Not incl in Mollys #####

#Calculate mean depth per individual
vcftools --bcf ${varcallbcf} --depth --out ${outfile}

#Calculate site depth
#####how is this different from one below? this incl in molly's. This is sum and sumsq depth.
#####below is mean and var depth.
vcftools --bcf ${varcallbcf} --site-depth --out ${outfile}

#Calculate site quality
vcftools --bcf ${varcallbcf} --site-quality --out ${outfile}
    ##### Not incl in Mollys #####
    
#Calculate proportion of missing data per individual
vcftools --bcf ${varcallbcf} --missing-indv --out ${outfile}

#Calculate proportion of missing data per site
vcftools --bcf ${varcallbcf} --missing-site --out ${outfile}
    
#Calculate heterozygosity
vcftools --bcf ${varcallbcf} --het --out ${outfile}

#To view these stats in R, see link above.
echo "Prefilter stats has finished."
