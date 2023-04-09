#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Preparing for impuation:
# Imputing tuturuatu_all bcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call bcf that has been filtered for 
#   4x coverage, 0.2 site missingness, 10 minGQ and no LD.

#Environment: impute

sppdir=~/data/tuturuatu_all_rm_a09/
## Edit to be run specific

# Making a directory to hold the imputation work.
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}filtered_input/
subsetdir=${impdir}vcf_subsets/

# Setting the filtered bcf to be working on and making a working file copy
cp ${sppdir}bcf/filter_strand_bias/impute_sb_files/Tuturuatu_VariantCalls_0.2site_missing_MinGQ10_0.6SP.vcf.gz* ${impdir}filtered_input/
filtervcf=${impdir}filtered_input/Tuturuatu_VariantCalls_0.2site_missing_MinGQ10_0.6SP.vcf.gz
##  Edit to be chosen filter specific

# Extracting TLR contigs
    # Extract the TLR contigs from the chosen filtered bcf
    echo "Extracting TLR contigs"
    bcftools view ${filtervcf} -r jcf7180002669510,jcf7180002696225,jcf7180002687310,jcf7180002686685,\
        jcf7180002693511,jcf7180002688524,jcf7180002688267,jcf7180002696332,jcf7180002694481,jcf7180002693589 \
        -O z -o ${subsetdir}tuturuatu_tlr_contigs.vcf.gz
        ## This could also be a bcf with -O b and .bcf. Not sure which format to use yet.

    # Index the new TLR contig vcf
    echo ""; echo "Indexing the TLR contig vcf"
    bcftools index ${subsetdir}tuturuatu_tlr_contigs.vcf.gz -f

# Setting reference (high coverage) population and study (low coverage) population
    echo ""; echo "Setting reference and study populations"
    # Extracting individual IDs
    bcftools query -l ${subsetdir}tuturuatu_tlr_contigs.vcf.gz > ${subsetdir}indv.ID

    # Subsetting study population
    grep -E "CR0[^6]|CR1|I16487" ${subsetdir}indv.ID > ${subsetdir}study.ID
        # This extracts indvs CR01-19, but not CR06. Also I16487.

    # Subsetting reference population
    grep -v -f ${subsetdir}study.ID ${subsetdir}indv.ID > ${subsetdir}ref.ID
        # This extracts indvs not present in the study population.

# Subsetting the tlr contig bcf into reference and study populations
    echo ""; echo "Subsetting the tlr contig bcf into reference and study populations"
    # Reference population
    echo "Subsetting reference population TLR contig file"
    bcftools view -O z -o ${subsetdir}tuturuatu_ref.vcf.gz -S ${subsetdir}ref.ID ${subsetdir}tuturuatu_tlr_contigs.vcf.gz
    # Index the ref TLR contig vcf
    echo "Indexing reference population TLR contig file"
    bcftools index ${subsetdir}tuturuatu_ref.vcf.gz -f

    # Study population
    echo ""; echo "Subsetting study population TLR contig file"
    bcftools view -O z -o ${subsetdir}tuturuatu_study.vcf.gz -S ${subsetdir}study.ID ${subsetdir}tuturuatu_tlr_contigs.vcf.gz
    # Index the study TLR contig vcf
    echo "Indexing study population TLR contig file"
    bcftools index ${subsetdir}tuturuatu_study.vcf.gz -f

wait

# Phasing 
## Is this necessary?
    # Phasing the reference population
    echo ""; echo "Phasing and indexing the reference population TLR contig file"
    beagle gt=${subsetdir}tuturuatu_ref.vcf.gz out=${subsetdir}tuturuatu_ref_phased
    bcftools index ${subsetdir}tuturuatu_ref_phased.vcf.gz -f

## Should I be phasing the study popl?
    # Phasing the study population
    echo ""; echo "Phasing and indexing the study population TLR contig file"
    beagle gt=${subsetdir}tuturuatu_study.vcf.gz out=${subsetdir}tuturuatu_study_phased
    bcftools index ${subsetdir}tuturuatu_study_phased.vcf.gz -f

echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found (phased) at ${subsetdir}"
echo "Now ready for imputing!!"