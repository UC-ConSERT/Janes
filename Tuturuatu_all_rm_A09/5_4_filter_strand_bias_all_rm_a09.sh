#!/bin/bash -e

#07 March 2023
#Olivia Janes, from Molly Magid
#Filtering the filtered bcfs, to remove any sites where strand bias adjusted phred score is <60
#Must be converted to a vcf first.

sppdir=~/data/tuturuatu_all_rm_a09/
## Edit to be run specific

filterdir=${sppdir}bcf/filter_trial/

# Making a directory to hold the strand bias filtered bcfs.
mkdir -p ${sppdir}bcf/filter_strand_bias/ ${sppdir}bcf/filter_trial/vcf/
vcfdir=${sppdir}bcf/filter_trial/vcf/
sbiasdir=${sppdir}bcf/filter_strand_bias/

# Convert no LD filter files to vcf.gz format
for bcf in ${filterdir}noLD/*bcf
do
	base=$(basename ${bcf} .bcf)
    # Convert to vcf.gz and index
	bcftools view ${bcf} --write-index -O z -o ${vcfdir}${base}.vcf.gz
	echo "finished ${base}"
	echo ""
done

# Convert LD filter files to vcf.gz format
for bcf in ${filterdir}LD_filter/*bcf
do
	base=$(basename ${bcf} .bcf)
    # Convert to vcf.gz and index
	bcftools view ${bcf} --write-index -O z -o ${vcfdir}${base}.vcf.gz
	echo "finished ${base}"
	echo ""
done

# Convert impute filter files to vcf.gz format
for bcf in ${filterdir}impute/*bcf
do
	base=$(basename ${bcf} .bcf)
    # Convert to vcf.gz and index
	bcftools view ${bcf} --write-index -O z -o ${vcfdir}${base}.vcf.gz
	echo "finished ${base}"
	echo ""
done

# Filter bcf files for strand bias.
    #This sets individual sites with SP <60 to "."

for file in ${vcfdir}*.vcf.gz
do
    echo "Filtering ${file} for SP <60"
    base=$(basename ${file} .vcf.gz)
    bcftools +setGT ${file} -- -t q -n . -i 'FORMAT/SP>60' > ${sbiasdir}${base}_0.6SP.vcf
done



echo "Filtering for strand bias is complete."

