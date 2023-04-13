#!/bin/bash -e

# 13 April 2023
# Olivia Janes
# Preparing for impuation:
# Imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters.

#Environment: impute

sppdir=~/data/tuturuatu_all_vcf/impute/validation/
    ## Edit to be run specific
beaglejar=~/data/programs/beagle.22Jul22.46e.jar
    ##Define location of beagle 5.4 program.
    ##Beagle can be downloaded using: wget http://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
study_list="A09|A11_nodup.bam|B10_nodup.bam|CR20_nodup.bam|CT07_nodup.bam|CT11_nodup.bam|E10_nodup.bam|\
                F09_nodup.bam|I16468_nodup.bam|I16476_nodup.bam"

# Making a directory to hold the imputation work.
mkdir -p ${sppdir}impute/
impdir=${sppdir}impute/
mkdir -p ${impdir}vcf_subsets/ ${impdir}vcf_finals/
subsetdir=${impdir}vcf_subsets/
finaldir=${impdir}vcf_finals/


# Extracting TLR contigs
    for vcf in ${sppdir}bcf/filter_trial/impute/*vcf.gz
    do
        base=$(basename ${vcf} _0.6SP.vcf.gz)
        # Extract the TLR contigs from the filtered vcf files
        # This is done in 3 files as bcftools only seems to be able to handle 4 "windows"/contigs at a time
        echo "Extracting TLR contigs"
        bcftools view --threads 16 ${vcf} -r jcf7180002669510,jcf7180002696225,jcf7180002687310,jcf7180002686685 \
            -O z -o ${subsetdir}${base}_1_tlr_contigs.vcf.gz
        bcftools view --threads 16 ${vcf} -r jcf7180002693511,jcf7180002688524,jcf7180002688267,jcf7180002696332 \
            -O z -o ${subsetdir}${base}_2_tlr_contigs.vcf.gz
        bcftools view --threads 16 ${vcf} -r jcf7180002694481,jcf7180002693589 \
            -O z -o ${subsetdir}${base}_3_tlr_contigs.vcf.gz
    done
    wait

# Index the new TLR contig vcf
    echo ""; echo "Indexing the TLR contig vcfs"
    for file in ${subsetdir}*_*_tlr_contigs.vcf.gz
    do
        bcftools index -f --threads 16 ${file}
    done

# Setting reference (high coverage) population and study (low coverage) population
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

# Phasing the reference panel
    for file in ${subsetdir}*_ref.vcf.gz
    do
        base=$(basename ${file} .vcf.gz)
        echo ""; echo "Phasing and indexing ${base} TLR contig file"
        java -jar ${beaglejar} gt=${file} out=${finaldir}${base}_phased
        bcftools index -f --threads 16 ${finaldir}${base}_phased.vcf.gz
    done


echo ""
echo "Script has finished preparing filtered variant call file ${filtervcf}."
echo "Your final study and ref files can be found at ${finaldir}"
echo "Now ready for imputing!!"