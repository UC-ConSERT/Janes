#!/bin/bash -e 
set -e

#14 April 2023

#Olivia Janes adapted from Molly Magid and Jana Wold
#Tuturuatu variant calling and preparing files for filtering
#Using downsampled files for imputation validation

##### Must be edited to be sample specific #####
sppdir=~/data/tuturuatu_all_vcf/
nodupbamdir=~/data/tuturuatu_all/nodup_bam/
        #Location of non-downsampled bam files to variant call alongside the downsampled bams
ref=~/data/tuturuatu_all/ref_genome/Maui_merged_assembly.fa
         #reference genome for alignment
#Ensure that the tlr_contigs.bed is in the ds_bams_tlr/ directory
        
#Preparing the directories
        dsdir=${sppdir}impute/downsampling_trial/
        mkdir -p ${nodupbamdir}not_var_calling ${dsdir}ds_bams_tlr/
        dsbamdir=${dsdir}ds_bams/
                #directory that holds the downsampled bams.
        dsbamtlrdir=${dsdir}ds_bams_tlr/
                #directory that holds the downsampled bams, only tlr contig regions.
        scriptdir=~/data/general_scripts/
        species="Tuturuatu"

#Setting up directories for each downsample trial
        mkdir -p ${dsdir}0.05_downsample/ ${dsdir}0.1_downsample/ ${dsdir}0.2_downsample/ ${dsdir}0.3_downsample/

        for ds in {0.05,0.1,0.2,0.3}
        do
                folder=${dsdir}${ds}_downsample/
                mkdir -p ${folder}bcf_tlr_${ds}/ ${folder}chunks_tlr_${ds}/
        done



#Remove the downsample original bams from the nodup folder to ensure they are not variant called
        # Define the list of files to downsample
        file_list="A09_nodup.bam A11_nodup.bam B10_nodup.bam CR20_nodup.bam CT07_nodup.bam CT11_nodup.bam E10_nodup.bam \
                F09_nodup.bam I16468_nodup.bam I16476_nodup.bam"
        echo "Moving files not to be variant called"
        for file in ${file_list}
        do
                mv ${nodupbamdir}${file}* ${nodupbamdir}not_var_calling/ || true
                        #|| true : this means that if this doesn't find a file to move, it won't stop the whole script
        done

#Subset the bam files so that they only contain the TLR contigs
        #Extracting TLR region for downsampled bam
        for file in ${dsbamdir}*.bam
        do
                echo "Extracting TLR region for ${file}"
                base=$(basename ${file} .bam)
                samtools view -@ 16 -b -h ${file} -o ${dsbamtlrdir}${base}_tlrs.bam -M -L ${dsbamtlrdir}tlr_contigs.bed
        done
        
        #Extracting TLR region for all other samples
        for file in ${nodupbamdir}*.bam
        do
                echo "Extracting TLR region for ${file}"
                base=$(basename ${file} .bam)
                samtools view -@ 16 -b -h ${file} -o ${dsbamtlrdir}${base}_tlrs.bam -M -L ${dsbamtlrdir}tlr_contigs.bed
        done

#Variant calling for each downsampling trial.
for ds in {0.05,0.1,0.2,0.3}
do
        echo""; echo ""; echo "Variant calling for ${ds} downsampling trial"
        #chunk bam files for mpileup
                echo ""; echo "Chunking ${ds} downsampling files for mpileup"
                ls ${dsbamtlrdir}*_nodup_tlrs.bam > ${dsbamtlrdir}${species}_${ds}ds_bam_list.txt
                ls ${dsbamtlrdir}*_${ds}subset_downsampled_tlrs.bam >> ${dsbamtlrdir}${species}_${ds}ds_bam_list.txt

                chunksdir=${dsdir}${ds}_downsample/chunks_tlr_${ds}/
                perl ${scriptdir}split_bamfiles_tasks.pl \
                        -b ${dsbamtlrdir}${species}_${ds}ds_bam_list.txt \
                        -g $ref -n 16 -o ${chunksdir} | parallel -j 16 {}

        #run mpileup on chunks of bam files
                echo ""; echo "Running mpileup on chunks of bam files"
                bcf_file=${dsdir}${ds}_downsample/bcf_tlr_${ds}/
                        #bcf file output

                for ((i=1; i<=16; i++))
                do
                        bcftools mpileup \
                                --threads 16 \
                                -f $ref \
                                -a AD,ADF,ADR,DP,SP \
                                -O b -o ${bcf_file}${species}_${i}_raw.bcf \
                                ${chunksdir}${i}/* &
                done
                wait
                echo "mpileup is done running." 


        #variant calling on bcf files
                echo ""; echo "Beginning variant calling..."
                for file in ${bcf_file}*.bcf
                do
                base=$(basename $file .bcf)
                bcftools call --threads 16 $file -mv -O b -f GQ -o ${bcf_file}${base}_tlr_VariantCalls.bcf &    
                done
                wait
                echo "Variant calling is complete." 


        #prepare files for filtering
                echo ""; echo "Preparing files for filtering..."
                for file in ${bcf_file}*Calls.bcf
                do
                        base=$(basename $file _VariantCalls.bcf)
                        #reheader each chunked bcf so it has the same sample names
                        bcftools reheader -s ${dsbamtlrdir}${species}_${ds}ds_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf  
                        wait
                        #put bcf files names into a list for concatenation
                        ls ${bcf_file}${base}_reheader.bcf >> ${bcf_file}list_of_bcf.txt 
                done
                echo "Preparing files complete. Concatenating chunked bcf files"

        #concatenate the chunked bcf files
                bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O b -o ${bcf_file}${species}_VariantCalls_${ds}ds_concat.bcf --threads 16
                echo ""; echo "bcf file is ready for filtering!"

done

echo ""; echo ""; echo "Entire script is complete!"