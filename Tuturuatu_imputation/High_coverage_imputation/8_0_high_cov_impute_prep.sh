#!/bin/bash -e
set -e

# 13 April 2023
# Olivia Janes
# Preparing for impuation:
# Imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters.

#Environment: impute

mkdir -p ~/data/tuturuatu_all_vcf/impute/truth/
sppdir=~/data/tuturuatu_all_vcf/impute/truth/
    ## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Making a directory to hold the imputation work.
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}vcf_finals/
subsetdir=${impdir}vcf_subsets/
finaldir=${impdir}vcf_finals/

study_list="A09|A11|B10|CR20|CT07|CT11|E10|F09|I16468|I16476"


#Copying over the TLR contig subsetted vcfs from the low coverage trial
    cp ~/data/tuturuatu_all_vcf/impute/vcf_subsets/*_*_tlr_contigs.vcf* ${subsetdir}

# Setting reference (high coverage) population and study (high coverage validation) population
    echo ""; echo "Setting reference and study populations"
    for file in ${subsetdir}*_1_tlr_contigs.vcf.gz
    do
        # Extracting individual IDs
        bcftools query -l ${file} > ${subsetdir}indv.ID

        # Subsetting removal (old study) population
        grep -E "CR0[^6]|CR1|I16487" ${subsetdir}indv.ID > ${subsetdir}remove.ID
            # This extracts indvs CR01-19, but not CR06. Also I16487.
        
        # Subsetting study population
        grep -E "${study_list}" ${subsetdir}indv.ID > ${subsetdir}study.ID
        grep -E "${study_list}" ${subsetdir}indv.ID >> ${subsetdir}remove.ID

        # Subsetting reference population
        grep -v -f ${subsetdir}remove.ID ${subsetdir}indv.ID > ${subsetdir}ref.ID
            # This extracts indvs not present in the study population.
        break
            # This stops the loop after one iteration.
    done


# Subsetting the tlr contig bcf into reference and study populations
    for file in ${subsetdir}*_*_tlr_contigs.vcf.gz
    do
        base=$(basename ${file} _tlr_contigs.vcf.gz)
        echo ""; echo "Subsetting ${base} TLR contig file into ref and study popls."
        # Reference population
        bcftools view -O z -o ${subsetdir}${base}_ref.vcf.gz -S ${subsetdir}ref.ID ${file}
        # Index the ref TLR contig vcf
        bcftools index ${subsetdir}${base}_ref.vcf.gz -f --threads 16

        # Study population
        bcftools view -O z -o ${finaldir}${base}_study.vcf.gz -S ${subsetdir}study.ID ${file}
        # Index the study TLR contig vcf
        bcftools index ${finaldir}${base}_study.vcf.gz -f --threads 16
    done


# Phasing the reference panel, trialling different Ne (as beagle imputes the reference panel during phasing)
    for file in ${subsetdir}*_ref.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Phasing and indexing ${base} TLR contig file"

        for test_ne in {50,100,500}
        do
            java -jar ${beaglejar} gt=${file} ne=${test_ne} em=false nthreads=16 \
                out=${finaldir}${base}_${test_ne}ne_phased 
            bcftools index -f --threads 16 ${finaldir}${base}_${test_ne}ne_phased.vcf.gz
        done
    done

# Phasing WITHOUT setting Ne, for comparison
    for file in ${subsetdir}*_ref.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Phasing and indexing ${base} TLR contig file"

        java -jar ${beaglejar} gt=${file} em=true nthreads=16 \
            out=${finaldir}${base}_defaultne_phased 
        bcftools index -f --threads 16 ${finaldir}${base}_defaultne_phased.vcf.gz
    done



echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
echo "Now ready for imputing!!"