## Assessing the Potential for Genetic Rescue to Improve Disease Outcomes in One of the World's Rarest Shorebirds

# Chapter 2: Characterising immune gene diversity in translocated wild tchūriwat’ (*Thinornis novaeseelandiae*)


### Programmes used
#### Environment: 
fastqc
- f
samtools

### 1_Prepping raw files:
I aligned these trimmed fastq files to the reference genome using the Burrows-Wheeler 
Alignment tool v0.7.17 to produce whole genome alignments (Li & Durbin, 2009). I then used 
SAMtools v1.17 to convert the BAM files to SAM files and then sort and index them (Li et al., 
2009). SAM files were tidied by filling in mate coordinates and inserting size fields (fixmate) and 
then marking and removing duplicates (markdup) with SAMtools (Li et al., 2009). Qualimap 
v2.2.2, Mosdepth v0.3.3 and SAMtools (flagstat) were used to assess the quality of samples 
before and after the tidy-up and ensure that duplicates had been removed (García-Alcalde et al., 
2012; Li et al., 2009; Pedersen & Quinlan, 2018)