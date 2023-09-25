## Assessing the Potential for Genetic Rescue to Improve Disease Outcomes in One of the World's Rarest Shorebirds

# Chapter 2: Characterising immune gene diversity in translocated wild tchūriwat’ (*Thinornis novaeseelandiae*)


### Programmes used:
#### Environment: 
fastqc
- f

samtools

bcftools

### 1_preparing_raw_files:
As part of this aligned project, sequences were quality checked and
trimmed (see Magid et al. (2022) for methods used). 


### 2_alignment:
I aligned these trimmed fastq files to the reference genome using the Burrows-Wheeler 
Alignment tool v0.7.17 to produce whole genome alignments (Li & Durbin, 2009). I then used 
SAMtools v1.17 to convert the BAM files to SAM files and then sort and index them (Li et al., 
2009). 


### 3_alignment_tidyup:
SAM files were tidied by filling in mate coordinates and inserting size fields (fixmate) and 
then marking and removing duplicates (markdup) with SAMtools (Li et al., 2009). Qualimap 
v2.2.2, Mosdepth v0.3.3 and SAMtools (flagstat) were used to assess the quality of samples 
before and after the tidy-up and ensure that duplicates had been removed (García-Alcalde et al., 
2012; Li et al., 2009; Pedersen & Quinlan, 2018).


### 4_variant_calling:


### 5_filtering:

### 6_imputation:
1_Imputation:
- Imputation of the low coverage (<5x) samples.


#### Trial names
[Not important, just for record keeping.]

Starts with tuturuatu, tuturuatu_trial_2 & tuturuatu_all (tuturuatu_all is the most recent, but some of the other scripts are needed for initial fq.gz & bam file prep).

Tuturuatu_all_rm_A09 & tuturuatu_all_rm_bad were trials  that figured out there were issues with the bcf file format downstream from variant calling, and tried variant calling with subsets of the data to see if there was an improvement in quality (in summary: ignore these scripts).  

Tuturuatu_all_vcf is the updated version of tuturuatu_all where issues with the vcf/bcf file format were resolved, starting with 5_0_filtering. This continues with post-filtering stats, where it was realised that the 2022 LIC files were too low quality to use. Then, the imputation route began under trial tuturuatu_imputation, starting again with 5_0_preimputation_filtering.

#### References: