#!/bin/bash -e

# 04 April 2023
# Olivia Janes
# Preparing for impuation:
# Imputing tuturuatu_all bcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call bcf that has been filtered for 
#   4x coverage, 0.2 site missingness, 10 minGQ and no LD.

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/
chosenfilter=Tuturuatu_VariantCalls_4x_coverage_0.2site_missing_MinGQ10_0.6SP
## Edit to be run specific

# Making a directory to hold the imputation work.
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}filtered_input/ ${impdir}vcf_finals/
subsetdir=${impdir}vcf_subsets/
finaldir=${impdir}vcf_finals/

# Setting the filtered bcf to be working on and making a working file copy
#cp ${sppdir}bcf/filter_strand_bias/${chosenfilter}.vcf.gz* ${impdir}filtered_input/
filtervcf=${impdir}filtered_input/${chosenfilter}.vcf.gz
##  Edit to be chosen filter specific

# Extracting TLR contigs
    # Extract the TLR contigs from the chosen filtered bcf
    # This is done in 3 files as bcftools only seems to be able to handle 4 "windows"/contigs at a time
    echo "Extracting TLR contigs"
    bcftools view --threads 16 ${filtervcf} -r jcf7180002669510,jcf7180002696225,jcf7180002687310,jcf7180002686685, \
        -O z -o ${subsetdir}tuturuatu_1_tlr_contigs.vcf.gz
    bcftools view --threads 16 ${filtervcf} -r jcf7180002693511,jcf7180002688524,jcf7180002688267,jcf7180002696332 \
        -O z -o ${subsetdir}tuturuatu_2_tlr_contigs.vcf.gz
    bcftools view --threads 16 ${filtervcf} -r jcf7180002694481,jcf7180002693589 \
        -O z -o ${subsetdir}tuturuatu_3_tlr_contigs.vcf.gz

    # Index the new TLR contig vcf
    echo ""; echo "Indexing the TLR contig vcfs"
    for file in ${subsetdir}tuturuatu_*_tlr_contigs*
    do
        bcftools index -f --threads 16 ${file}
    done

# Setting reference (high coverage) population and study (low coverage) population
    echo ""; echo "Setting reference and study populations"
    # Extracting individual IDs
    bcftools query -l ${subsetdir}tuturuatu_1_tlr_contigs.vcf.gz > ${subsetdir}indv.ID

    # Subsetting study population
    grep -E "CR0[^6]|CR1|I16487" ${subsetdir}indv.ID > ${subsetdir}study.ID
        # This extracts indvs CR01-19, but not CR06. Also I16487.

    # Subsetting reference population
    grep -v -f ${subsetdir}study.ID ${subsetdir}indv.ID > ${subsetdir}ref.ID
        # This extracts indvs not present in the study population.


# Phasing    
    for file in ${subsetdir}tuturuatu_*_tlr_contigs.vcf.gz
    do
        base=$(basename ${file} _tlr_contigs.vcf.gz)
        echo ""; echo "Phasing and indexing ${base} TLR contig file"
        beagle gt=${file} out=${subsetdir}${base}_phased
        bcftools index -f --threads 16 ${subsetdir}${base}_phased.vcf.gz
        wait

        # Subsetting the tlr contig bcf into reference and study populations
            echo ""; echo "Subsetting the ${base} tlr contig file into reference and study populations"
            # Reference population
            echo "Subsetting reference population TLR contig file"
            bcftools view --threads 16 -O z -o ${finaldir}${base}_ref.vcf.gz -S ${subsetdir}ref.ID ${subsetdir}${base}_phased.vcf.gz
            wait
            # Index the ref TLR contig vcf
            echo "Indexing reference population TLR contig file"
            bcftools index --threads 16 ${finaldir}${base}_ref.vcf.gz

            # Study population
            echo ""; echo "Subsetting study population TLR contig file"
            bcftools view --threads 16 -O z -o ${finaldir}${base}_study.vcf.gz -S ${subsetdir}study.ID ${subsetdir}${base}_phased.vcf.gz
            wait
            # Index the study TLR contig vcf
            echo "Indexing study population TLR contig file"
            bcftools index -f --threads 16 ${finaldir}${base}_study.vcf.gz 
    done


echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
echo "Now ready for imputing!!"