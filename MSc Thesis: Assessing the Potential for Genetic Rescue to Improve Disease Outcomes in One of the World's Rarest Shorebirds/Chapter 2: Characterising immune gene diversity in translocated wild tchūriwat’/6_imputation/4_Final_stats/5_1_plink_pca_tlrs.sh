#!/bin/bash -e 
set -e

# 09 May 2023

# Olivia Janes
# Running a plink pca on the TLR contigs.
# From: tuturuatu_imputation

## Environment: plink

# Setting up
    ##  Needs to be edited to be run specific   ##
        sppdir=~/data/tuturuatu_all_vcf/
        tlrdir=${sppdir}final_outputs/final_vcfs/
        tlrvcf=${tlrdir}Tuturuatu_VariantCalls_0.9miss_low_cov_final.vcf.gz

    #Setting directories
        mkdir -p ${sppdir}final_outputs/final_stats/pca_plink/tlr_pca/ ${tlrdir}LD_filtered_vcf/
        lddir=${tlrdir}LD_filtered_vcf/
        statsdir=${sppdir}final_outputs/final_stats/pca_plink/tlr_pca/


#Perform linkage filtering
    echo "Creating a list of sites to be LD filtered"
    plink --vcf ${tlrvcf} --double-id --allow-extra-chr \
        --set-missing-var-ids @:# \
        --indep-pairwise 50 10 0.1 --out ${lddir}LD_filtered_sites

#Create PCAs
    echo "Filtering and creating PCA"
    plink --vcf ${tlrvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --extract ${lddir}LD_filtered_sites.prune.in \
        --make-bed --pca --out ${statsdir}Tuturuatu_VariantCalls_tlrimp_LD_filtered

    echo ""; echo "Creating non-filtered PCA"
    plink --vcf ${tlrvcf} --double-id --allow-extra-chr --set-missing-var-ids @:# \
        --make-bed --pca --out ${statsdir}Tuturuatu_VariantCalls_tlrimp

echo ""; echo "Script is complete"
echo "Download PCA files with:"
echo "rsync -rav rccuser:/home/rccuser/data/tuturuatu_all_vcf/final_outputs/final_stats/pca_plink/tlr_pca ./"
    

