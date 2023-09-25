#!/bin/bash -e 
set -e

# 09 May 2023

# Olivia Janes
# Running a plink pca on the whole genome vcf created in 4_0_wgs_genotypes_for_pca.sh.
# Also running a plink pca for the original (og) wgs vcf, without the imputed TLR region.
# From: tuturuatu_imputation

## Environment: plink

# Setting up
    ##  Needs to be edited to be run specific   ##
        sppdir=~/data/tuturuatu_all_vcf/
        wgsdir=${sppdir}final_outputs/final_vcfs/whole_genome/
        wgsvcf=${wgsdir}Tuturuatu_VariantCalls_wgs_tlrimp_merged.vcf.gz
        ogvcf=${sppdir}bcf/filter_trial/noLD/Tuturuatu_VariantCalls_5x_coverage_0.2site_missing_MinGQ10.vcf.gz

    #Setting directories
        mkdir -p ${sppdir}final_outputs/final_stats/pca_plink/og_pca/ ${wgsdir}LD_filtered_vcf/
        lddir=${wgsdir}LD_filtered_vcf/
        statsdir=${sppdir}final_outputs/final_stats/pca_plink/


#Perform linkage filtering
    echo "Creating a list of sites to be LD filtered"
    plink --vcf ${wgsvcf} --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --indep-pairwise 50 10 0.1 --out ${lddir}LD_filtered_sites
    
    plink --vcf ${ogvcf} --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --indep-pairwise 50 10 0.1 --out ${lddir}LD_filtered_sites_og   

#Create PCAs
    echo "Filtering and creating PCA"
    plink --vcf ${wgsvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --extract ${lddir}LD_filtered_sites.prune.in \
        --make-bed --pca --out ${statsdir}Tuturuatu_VariantCalls_wgs_tlrimp_LD_filtered

    echo ""; echo "Creating non-filtered PCA"
    plink --vcf ${wgsvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${statsdir}Tuturuatu_VariantCalls_wgs_tlrimp

    echo ""; echo "Creating PCA based on non-imputed WGS vcf, filtered and not filtered for LD"
    plink --vcf ${ogvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --extract ${lddir}LD_filtered_sites_og.prune.in \
        --make-bed --pca --out ${statsdir}og_pca/Tuturuatu_VariantCalls_wgs_notlr_LD_filtered

    plink --vcf ${ogvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${statsdir}og_pca/Tuturuatu_VariantCalls_wgs_notlr   

echo ""; echo "Script is complete"
echo "Download PCA files with:"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/pca_plink ./"
    

