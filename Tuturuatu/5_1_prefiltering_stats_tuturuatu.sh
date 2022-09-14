#!/bin/sh

#12 Aug 2022
#Adapted by Olivia Janes from https://speciationgenomics.github.io/filtering_vcfs/
#Pre filtering statistics calculated to guide vcf filtering and compare to post-filtering statistics.


sppdir=~/data/tuturuatu/
bcfdir=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

mkdir -p ${bcfdir}stats/ ${bcfdir}stats/prefilter_stats/

varcallbcf=${bcfdir}Tuturuatu_VariantCalls_concat.bcf
    #output concatenated file from variant calling
    ##### Must be edited to be sample specific #####
    ##### Must change to bcf or vcf as previous file was #####
outfile=${bcfdir}stats/prefilter_stats/Tuturuatu_VariantCalls_prefilter
    #defining the output file prefix

#If uncompressed vcf, should be compressed here using bgzip [variantcalls].vcf
#Index bcf
echo "Indexing variant calls file"
bcftools index ${varcallbcf}
wait

#Calculate stats
echo "Calculating pre-filter stats"
    ##### Check these stats are the same as calculated afterwards #####
vcftools --bcf ${varcallbcf} --freq2 --out ${outfile} --max-alleles 2
    ##### output may be changed to --gzvcf so check ######
    #Calculate allele frequency for each variant
vcftools --bcf ${varcallbcf} --depth --out ${outfile}
    #Calculate mean depth per individual
vcftools --bcf ${varcallbcf} --site-mean-depth --out ${outfile}
    #Calculate mean depth per site
vcftools --bcf ${varcallbcf} --site-quality --out ${outfile}
    #Calculate site quality
vcftools --bcf ${varcallbcf} --missing-indv --out ${outfile}
    #Calculate proportion of missing data per individual
vcftools --bcf ${varcallbcf} --missing-site --out ${outfile}
    #Calculate proportion of missing data per site
vcftools --bcf ${varcallbcf} --het --out ${outfile}
    #Calculate heterozygosity

#To view these stats in R, see link above.
echo "Prefilter stats has finished."
