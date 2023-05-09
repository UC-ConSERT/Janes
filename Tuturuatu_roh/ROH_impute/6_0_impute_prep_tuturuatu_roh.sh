#!/bin/bash -e

# 05 May 2023
# Olivia Janes
# Preparing for imputation:
# Imputing tuturuatu_all vcf to improve SNP calling in low coverage (<5x) individuals.
# Imputation will be on the variant call vcf that has been through various test filters (5x and 10x).

#Environment: samtools

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

<<"COMMENTS"
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
COMMENTS


# Filter reference panel for missingness, and output removed sites into a text file. 

        for file in ${subsetdir}*x_ref.vcf.gz
        do
            echo ""; echo "Printing sites with present in less than 0.8 of indv."
            base=$(basename ${file} _ref.vcf.gz)
            vcftools --gzvcf ${file} \
                --max-missing 0.8 \
                --out ${subsetdir}${base}_ref_0.8miss \
                --removed-sites
            echo "Finished printing missing sites, output found at: ${subsetdir}"
        done


# Remove the missing sites from the vcfs.
    for i in {5,10}
    do
        for file in ${subsetdir}*${i}x_ref.vcf.gz
        do
            base=$(basename ${file} _ref.vcf.gz)

            #Ref
            vcftools --gzvcf ${file} \
                --exclude-positions ${subsetdir}Tuturuatu_VariantCalls_${i}x_ref_0.8miss.removed.sites \
                --recode \
                --recode-INFO-all \
                --out ${subsetdir}${base}_0.8miss_ref.vcf

            #Study
            vcftools --gzvcf ${finaldir}${base}_study.vcf.gz \
                --exclude-positions ${subsetdir}Tuturuatu_VariantCalls_${i}x_ref_0.8miss.removed.sites \
                --recode \
                --recode-INFO-all \
                --out ${finaldir}${base}_0.8miss_study.vcf

            echo "Renaming filter file to remove '.recode.vcf'"
            mv -i ${subsetdir}${base}_0.8miss_ref.vcf.recode.vcf ${subsetdir}${base}_0.8miss_ref.vcf
            mv -i ${finaldir}${base}_0.8miss_study.vcf.recode.vcf ${finaldir}${base}_0.8miss_study.vcf
        done
    done


# Convert filter files to vcf.gz format and index
    for vcf in ${subsetdir}*0.8miss_ref.vcf
    do
        base=$(basename ${vcf} _ref.vcf)
        
        #Ref
        echo "Converting ${vcf} to vcf.gz format"
        bcftools view ${vcf} -O z -o ${subsetdir}${base}_ref.vcf.gz --threads 16

        echo "Indexing ${vcf}.gz"
        bcftools index ${subsetdir}${base}_ref.vcf.gz --threads 16

        if [ -e "${subsetdir}${base}_ref.vcf.gz" ]; then
            echo "Removing ${subsetdir}${base}_ref.vcf"
            rm "${subsetdir}${base}_ref.vcf"
        fi

        #Study
        echo "Converting ${base}_study to vcf.gz format"
        bcftools view ${finaldir}${base}_study.vcf -O z -o ${finaldir}${base}_study.vcf.gz --threads 16

        echo "Indexing ${base}_study.vcf"
        bcftools index ${finaldir}${base}_study.vcf.gz --threads 16

        if [ -e "${finaldir}${base}_study.vcf.gz" ]; then
            echo "Removing ${finaldir}${base}_study.vcf"
            rm "${finaldir}${base}_study.vcf"
        fi

        echo ""
    done


#Subsetting into contigs
counter=0
chromosomes=$(cat ${subsetdir}contigs.txt)
for i in {5,10}
do
    for chrom in ${chromosomes}
    do
        for vcf in ${subsetdir}Tuturuatu_VariantCalls_${i}x_0.8miss_ref.vcf.gz
        do
            base=$(basename ${vcf} _ref.vcf.gz)
            # Extract the contigs from the filtered vcf files
            counter=$[counter+1]

            echo "Extracting ${chrom} for ${base} ref, counter is ${counter}"
            bcftools view --threads 16 ${vcf} -r ${chrom} \
                -O z -o ${subsetdir}${base}_${counter}_ref.vcf.gz
            echo "Indexing ${chrom} for ${base} ref"
            bcftools index -f --threads 16 ${subsetdir}${base}_${counter}_ref.vcf.gz

            echo "Extracting ${chrom} for ${base} study, counter is ${counter}"
            bcftools view --threads 16 ${finaldir}${base}_study.vcf.gz -r ${chrom} \
                -O z -o ${finaldir}${base}_${counter}_study.vcf.gz
            echo "Indexing ${chrom} for ${base} study"
            bcftools index -f --threads 16 ${finaldir}${base}_${counter}_study.vcf.gz
        done
    done
done

# Phasing the reference panel
    for file in ${subsetdir}*x_0.8miss_*_ref.vcf.gz
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