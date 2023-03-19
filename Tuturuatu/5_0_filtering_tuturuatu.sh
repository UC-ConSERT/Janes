#!/bin/bash -e 

#31 Aug 2022
#Molly Magid adapted by Olivia Janes
#Filtering trials for Tuturuatu variant calls, from [Molly's github](https://github.com/UC-ConSERT/Magid_et_al/blob/main/3_filtering.sh)

sppdir=~/data/tuturuatu/
work=${sppdir}bcf/
    #directory where files to filter are ("bcf_file" from script 4_variant_calling_filter_prep_tuturuatu.sh)

mkdir -p ${work}filter_trial/ 
mkdir -p ${work}filter_trial/noLD/ ${work}filter_trial/LD_filter/ ${work}/stats/

vcf_out=${work}filter_trial/
noLD=${vcf_out}noLD/
LD=${vcf_out}LD_filter/

#for loop to filter file with different values for parameters including
#missingness, depth, and GQ
for bcf in ${work}*_concat.bcf
do
    base=$(basename ${bcf} _concat.bcf)
    for i in {4..5} #filtering files for 4x and 5x depth
    do
        echo "Filtering SNPs for ${base}...." 
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_noMinGQ.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_noMinGQ.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ10.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ10.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 10 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.1site_missing_MinGQ20.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.9 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all &
        vcftools --bcf ${bcf} \
            --out ${noLD}${base}_${i}x_coverage_0.2site_missing_MinGQ20.bcf \
            --minDP ${i} \
            --maxDP 200  \
            --max-missing 0.8 \
            --maf 0.05 \
            --minQ 20 \
            --minGQ 20 \
            --remove-indels \
            --remove-filtered-all \
            --recode-bcf \
            --recode-INFO-all
    done;
done

# Removing '.recode.bcf' from noLD file names.
echo "Renaming noLD files to remove '.recode.bcf'"
for bcf in ${noLD}*.bcf
do
    echo ""
        echo "Renaming ${bcf}"
    name=$(basename $bcf .recode.bcf)
    rename "s/${bcf}/${name}/g" ${bcf}*
done 


for sample in ${datadir}I*_L003_val_1.fq.gz
do 
        echo $sample
        base=$(basename $sample _val_1.fq.gz)
    #base=I164xx-L1_Sxxx_L003
        echo $base
        name=$(echo $base | sed 's/-L1_S[0-9][0-9][0-9]_L003//g')
        name1=$(echo $name | sed 's/-L1_S[0-9][0-9]_L003//g')
    #name=I16xx
        echo $name1
    #replaces ${base} with ${name1}, for all samples starting with base (therefore includes R2 with it)
        rename "s/${base}/${name1}/g" ${datadir}${base}*
    #should now be I16xx_val_x.fastq.gz
done



echo "Filtering for Linkage parameters..."
#for loop to filter previous filtered files for linkage
    ###### OJ note - changed -l to -m as -l was not recognised #######
for bcf in ${noLD}*.bcf
do
  base=$(basename ${bcf} .bcf)
    echo "Running light LD pruning at 0.8 for ${base}...."
   bcftools +prune \
        -m 0.8 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.8LD.bcf \
        ${bcf} &
    echo "Running moderate LD pruning at 0.6 for ${base}...."
    bcftools +prune \
        -m 0.6 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.6LD.bcf \
        ${bcf} &
    echo "Running strong LD pruning at 0.4 for ${base}...."
    bcftools +prune \
        -m 0.4 \
        -w 1000 \
        -O b \
        -o ${LD}${base}_0.4LD.bcf \
        ${bcf}
done



#calculating statistics for no linkage filtered files
for file in ${noLD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --site-depth &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --depth &
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-site &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-indv &
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --het
done

#calculating statistics for linkage filtered files
for file in ${LD}*.bcf
do
    base=$(basename ${file} .bcf)
    echo "Calculating depth for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --site-depth &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --depth &
    echo "Calculating missingness for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-site &
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --missing-indv &
    echo "Calculating individual heterozygosity for ${base}..."
    vcftools --bcf ${file} \
        --out ${work}stats/${base} \
        --het
done

echo "Script has finished running. Now whether it worked or not is another question..."