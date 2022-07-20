#Prepare beginning.




#chunk bam files for mpileup
ls ${processedbamdir}*.aligned.sorted.bam > ${processedbamdir}${species}bam_list.txt
perl ~/data/general_scripts/split_bamfiles_tasks.pl -b ${processedbamdir}${species}bam_list.txt -g $ref -n 13 -o ${chunksdir} | parallel -j 12 {}

#run mpileup on chunks of bam files
for ((i=1; i<=12; i++)); do
        bcftools mpileup -O b -f $ref -a AD,ADF,ADR,DP,SP -o ${bcf_file}${species}_${i}_raw.bcf  ${chunksdir}${i}/* &
done
wait
echo “mpileup is done running”

COMMENTS

#variant calling on bcf files
#for file in ${bcf_file}*.bcf
#do
#base=$(basename $file .bcf)
#bcftools call $file -mv -O b -f GQ -o ${bcf_file}${base}_VariantCalls.bcf &
#done
#wait
#echo “variant calling is complete”

#Molly doesn't have either of the below 2 lines in her script. it appears the list of bcf is created below on line 131
#>${bcf_file}list_of_bcf.txt   #do i need this part? I think this creates the list_of_bcf.txt file
#ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt #And this writes all of the variant calls names to it. However,
        #this may be done below after reheader. As in Molly's script. Final verdict: I THINK remove this part.


#prepare files for filtering with bgzip and indexing
#for file in ${bcf_file}*.vcf
#do
#base=$(basename $file _raw_VariantCalls.vcf)
#bcftools reheader -s ${bamdir}${species}_bam_list.txt ${file} -o ${bcf_file}${base}_reheader.bcf
#wait
#>list_of_bcf.txt
#ls ${bcf_file}*_VariantCalls.bcf >> ${bcf_file}list_of_bcf.txt
#done

#concatenate the chunked vcf files
#bcftools concat --file-list ${bcf_file}list_of_bcf.txt -O v -o ${bcf_file}${species}_VariantCalls_concat.vcf --threads 16
#echo “bcf file is ready for filtering!”