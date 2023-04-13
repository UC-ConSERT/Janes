#!/bin/bash -e 

#10 April 2023
#Olivia Janes
#Compiling SNP counts for the TLR contigs from filtering to compare between filtering methods.
#Repeated for LD filtered, noLD filtered, strand bias (SP) filtered and imputation filtered files.
#The outputs TLR_contigs_SNP_counts.txt is set up to be easily imported into a spreadsheet (comma delimited).

## NOTE this script doesn't work because the bed file doesn't have chromosome (contig) start and end coordinates.
##  It also doesn't work putting the contig id into the bcftools line (maxes out at 4)

sppdir=~/data/tuturuatu_all_vcf/
bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
sbiasdir=${sppdir}bcf/filter_strand_bias/
#Define location of tlr_contigs.bed file in script. Should be in bcf/


## TLR contigs SNP counts ##
    echo "Counting TLR contig SNPs"
    echo "TLR contig SNPs,Number" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt

    #Prefilter SNP counts
        echo "Prefilter TLR contig SNP counts"
        cd ${bcfdir}

        #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
        for file in ${bcfdir}*VariantCalls_concat.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -r ${bcfdir}tlr_contigs.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt
        done

    # Filtered files should be indexed in 5_0_filtering.sh
    # LD filtered SNP counts
        echo "LD filtered TLR SNP counts"
        cd ${filterdir}LD_filter/

        for file in ${filterdir}LD_filter/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_contigs.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt
        done


    # no LD filtered SNP counts
        echo "no LD filtered TLR SNP counts"
        cd ${filterdir}noLD/

        for file in ${filterdir}noLD/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_contigs.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt
        done

    # Imputation filtered SNP counts
        echo "Imputation filtered TLR SNP counts"
        cd ${filterdir}impute/filters/

        for file in ${filterdir}impute/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_contigs.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt
        done

    # Strand bias filtered SNP counts
        echo "Strand bias filtered TLR SNP counts"
        cd ${sbiasdir}

        for file in ${sbiasdir}*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_contigs.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_contigs_SNP_counts.txt
        done

echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/TLR_contigs_SNP_counts.txt"