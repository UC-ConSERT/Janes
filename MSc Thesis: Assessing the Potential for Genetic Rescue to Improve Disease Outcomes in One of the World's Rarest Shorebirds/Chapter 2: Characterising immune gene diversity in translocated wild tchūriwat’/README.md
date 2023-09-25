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
A custom Perl script was used to split the BAM files into chunks to speed up the variant 
calling process (Moraga, 2018). BCFtools v1.17 was used to produce BCF files with annotations 
‘AD,ADF,ADR,DP,SP’ (mpileup), and then for variant calling (Danecek et al., 2021; Li, 2011). 
Chunked BCF files were concatenated back together with BCFtools (Li, 2011). 


### 5_filtering:
#### 5_filtering_no_impute:
I then converted BCF files into VCF files using BCFtools (Li, 2011). I conducted filtering trials on the VCF file using 
VCFtools v0.1.16 to determine the final filter parameters (Danecek et al., 2011). The optimal 
filtering parameters for the whole genome resequencing VCF were determined based on mean
depth and missingness at the site and individual level, as well as site quality and resulting 
number of SNPs. 

A filtering trial was conducted here.
The filtering was conducted with VCFtools and the final parameters used were: 
Phred score (variant quality) > 20, genotype quality (GQ) > 10, minimum depth > 5, maximum 
depth < 50, maximum missingness < 0.2 and minor allele frequency > 0.05. I then used BCFtools 
to filter for strand-bias adjusted Phred-score < 60. This produced a final whole genome VCF
(variant call file) where called SNPs remaining were of adequate quality for any subsequent 
whole genome analyses (for the use of this whole-genome VCF, see Chapter 3).

After extracting individual depth and missingness stats for the whole genome and for 
within the TLR region (TLR regions were identified in an aligned project) using VCFtools, it was 
apparent that 16 of the 18 translocated wild individuals and 1 captive individual had very poor 
coverage and low depth, both genome-wide and within the TLRs (Danecek et al., 2011; Magid, 
2021). These individuals had an average depth of 3.2-5.3x before filtering, and after filtering
with the parameters set above, were missing 40-76% of SNP sites (compared to high-coverage
individuals that had an average depth of 6-19x and ~5% of sites missing). Importantly, within 
the TLR regions, these low-coverage individuals were missing 14-71% of the TLR SNP sites
found in high-coverage individuals after minimum quality filtering. I did not lessen the filtering 
parameters set above, as this may introduce high rates of error. As the intention of this process 
was to characterise TLR diversity in the 18 translocated wild birds, imputation was necessary to
augment the TLR regions of the low-coverage individuals.


#### 5_filtering_impute:
I used different filtering parameters for imputation than were used for the whole-genome
VCF. The filtering parameters were: Phred score (variant quality) > 20, genotype quality (GQ) > 
10, maximum depth < 50, and bi-allelic sites only, using VCFtools (Danecek et al., 2011). To test 
the effect of a minimum depth filter I trialled filtering for a minimum depth of 4x and 5x with 
VCFtools (Danecek et al., 2011). I then used BCFtools to filter for strand-bias adjusted Phred-score < 60. Minor allele frequency (MAF) and missingness in sites were not filtered pre-imputation, to retain as many sites as possible to increase the accuracy of imputation (Hui et al., 
2020).


### 6_imputation:
#### 1_Imputation: 
Imputation of the low coverage (<5x) samples.

#### 2_Validation_and_truth_imputation: 
Downsampling a subset of high coverage samples to validate the accuracy of imputation.

#### 3_Downsampling_diff_cov_trial: 
Downsampling the high coverage subset to lower depths and higher missingness to trial the limits of imputation on low coverage samples.

#### 4_Final_stats: 
Final stats for the low coverage imputation samples, including extracting genotypes, preparing PCAs, etc...

____________________________________________________________________________________________________


#### Trial names
[Not important, just for record keeping.]

Starts with tuturuatu, tuturuatu_trial_2 & tuturuatu_all (tuturuatu_all is the most recent, but some of the other scripts are needed for initial fq.gz & bam file prep).

Tuturuatu_all_rm_A09 & tuturuatu_all_rm_bad were trials  that figured out there were issues with the bcf file format downstream from variant calling, and tried variant calling with subsets of the data to see if there was an improvement in quality (in summary: ignore these scripts).  

Tuturuatu_all_vcf is the updated version of tuturuatu_all where issues with the vcf/bcf file format were resolved, starting with 5_0_filtering. This continues with post-filtering stats, where it was realised that the 2022 LIC files were too low quality to use. Then, the imputation route began under trial tuturuatu_imputation, starting again with 5_0_preimputation_filtering.

#### References: