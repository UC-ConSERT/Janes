#!/bin/bash -e 

#10 April 2023
#Olivia Janes
#Compiling total SNP counts and TLR SNP counts from filtering to compare between filtering methods.
#Repeated for LD filtered, noLD filtered, strand bias (SP) filtered and imputation filtered files.
#The outputs SNP_counts.txt and TLR_SNP_counts.txt are set up to be easily imported into a spreadsheet (comma delimited).

sppdir=~/data/tuturuatu_all_vcf/
bcfdir=${sppdir}bcf/
filterdir=${sppdir}bcf/filter_trial/
sbiasdir=${sppdir}bcf/filter_strand_bias/
#Define location of tlr_regions.bed file in script. Should be in bcf/

## Total SNP counts ##
    echo "Counting TOTAL SNPs"
    echo "TOTAL SNPs,Number" >> ${bcfdir}stats/SNP_counts.txt
    # Prefiltered SNP counts
    #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script.
        echo "Prefiltered SNP counts"
        cd ${bcfdir}

        for file in ${bcfdir}*VariantCalls_concat.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/SNP_counts.txt
        done

    # LD filtered SNP counts
        echo "LD filtered SNP counts"
        cd ${filterdir}LD_filter/

        for file in ${filterdir}LD_filter/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/SNP_counts.txt
        done

    # no LD filtered SNP counts
        echo "no LD filtered SNP counts"
        cd ${filterdir}noLD/

        for file in ${filterdir}noLD/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/SNP_counts.txt
        done

    # Imputation filtered SNP counts
        echo "Imputation filtered SNP counts"
        cd ${filterdir}impute/

        for file in ${filterdir}impute/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/SNP_counts.txt
        done

    # Strand bias filtered SNP counts
        echo "Strand bias filtered SNP counts"
        cd ${sbiasdir}

        for file in ${sbiasdir}*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/SNP_counts.txt
        done


## TLR SNP counts ##
    echo "Counting TLR SNPs"
    echo "TLR SNPs,Number" >> ${bcfdir}stats/TLR_SNP_counts.txt

    #Prefilter SNP counts
        echo "Prefilter TLR SNP counts"
        cd ${bcfdir}

        #For bcf file that was indexed in previous script (5_1_prefiltering_stats), and comes in compressed format from 4_variant_calling script. 
        for file in ${bcfdir}*VariantCalls_concat.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_SNP_counts.txt
        done

    # Filtered files should be indexed in 5_0_filtering.sh
    # LD filtered SNP counts
        echo "LD filtered TLR SNP counts"
        cd ${filterdir}LD_filter/

        for file in ${filterdir}LD_filter/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_SNP_counts.txt
        done


    # no LD filtered SNP counts
        echo "no LD filtered TLR SNP counts"
        cd ${filterdir}noLD/

        for file in ${filterdir}noLD/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_SNP_counts.txt
        done

    # Imputation filtered SNP counts
        echo "Imputation filtered TLR SNP counts"
        cd ${filterdir}impute/filters/

        for file in ${filterdir}impute/*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_SNP_counts.txt
        done

    # Strand bias filtered SNP counts
        echo "Strand bias filtered TLR SNP counts"
        cd ${sbiasdir}

        for file in ${sbiasdir}*.vcf.gz
        do
            base=$(basename ${file} .vcf.gz)
            snp=$(bcftools query -R ${bcfdir}tlr_regions.bed -f '%POS\n' ${file} | wc -l) 
            echo "${base},${snp}" >> ${bcfdir}stats/TLR_SNP_counts.txt
        done

echo "SNP counting is complete. Yay! Find output file at ${bcfdir}stats/SNP_counts.txt"