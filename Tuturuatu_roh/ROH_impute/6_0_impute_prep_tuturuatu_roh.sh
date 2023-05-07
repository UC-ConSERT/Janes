#!/bin/bash -e

# 05 May 2023
# Olivia Janes
# Preparing for imputation:
# Imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters (5x and 10x).

#Environment: impute

sppdir=~/data/tuturuatu_roh/
    ## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Setting the directories & making a directory to hold the imputation work.
filterdir=${sppdir}bcf/filter_trial/impute/
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}vcf_finals/
subsetdir=${impdir}vcf_subsets/
finaldir=${impdir}vcf_finals/


# Setting reference (high coverage) population and study (low coverage) population
    echo ""; echo "Setting reference and study populations"
    for file in ${filterdir}*.vcf.gz
    do
        # Extracting individual IDs
        bcftools query -l ${file} > ${subsetdir}indv.ID

        # Subsetting study population
        grep -E "CR0[^6]|CR1|I16487" ${subsetdir}indv.ID > ${subsetdir}study.ID
            # This extracts indvs CR01-19, but not CR06. Also I16487.

        # Subsetting reference population
        grep -v -f ${subsetdir}study.ID ${subsetdir}indv.ID > ${subsetdir}ref.ID
            # This extracts indvs not present in the study population.
        break
            # This stops the loop after one iteration.
    done


# Subsetting the vcf into reference and study populations
    for file in ${filterdir}*.vcf.gz
    do
        base=$(basename ${file} _0.6SP.vcf.gz)
        echo ""; echo "Subsetting ${base} file into ref and study popls."
        # Reference population
        bcftools view -O z -o ${subsetdir}${base}_ref.vcf.gz -S ${subsetdir}ref.ID ${file}
        # Index the ref contig vcf
        bcftools index ${subsetdir}${base}_ref.vcf.gz -f --threads 16

        # Study population
        bcftools view -O z -o ${finaldir}${base}_study.vcf.gz -S ${subsetdir}study.ID ${file}
        # Index the study contig vcf
        bcftools index ${finaldir}${base}_study.vcf.gz -f --threads 16
    done

# Phasing the reference panel
    for file in ${subsetdir}*_ref.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Phasing and indexing ${base} file"
        java -jar ${beaglejar} gt=${file} out=${finaldir}${base}_phased ne=100 em=false
        bcftools index -f --threads 16 ${finaldir}${base}_phased.vcf.gz
    done


echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
echo "Now ready for imputing!!"