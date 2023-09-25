# Assessing the Potential for Genetic Rescue to Improve Disease Outcomes in One of the World's Rarest Shorebirds
#### Masters Thesis submitted by Olivia Rose Janes
#### for Master of Science in Biology at the University of Canterbury
#### May 2023


Presented here are the scripts used for data analyses for the data chapters of my thesis. 
Provided below is a summary of the analyses conducted, along with a list of conda environments and relevant programs used.


## Chapter 2: Characterising immune gene diversity in translocated wild tchūriwat’ (*Thinornis novaeseelandiae*)
In this chapter, I charaterised toll-like receptor (TLR) single-nucleotide polymorphisms (SNPs) in captive and wild tchūriwat. I develop and test methods to impute low-coverage whole genome sequencing data when a large reference data set is not available.


## Chapter 3: Investigating associations between immune gene diversity and immune response in tchūriwat’



### Programs used:
#### Programs:
FastQC 
- Version 0.11.9

Trim-galore 
- Version 0.6.7 

Bcftools 
- Version 1.17

Rename 
- Version 1.601

PLINK 
- Version 1.90b6.21

Beagle
- Version 5.4

DnaSP 
- Version 6.12

PopART 
- Version 1.7

GATK 
- Version 4.4.0.0

Java 
- Version 17.0.6




#### Environments: 
fastqc
- fastqc
- trim-galore

mosdepth
- mosdepth
- samtools
- qualimap

samtools
- samtools
- bcftools
- vcftools
- rename

bcftools
- bcftools
- rename

impute
- beagle v5.2 (if using beagle v5.4, use: `java -jar ~/data/programs/beagle.22Jul22.46e.jar`)
- bcftools


[Link text](link address)

# Heading 1
## Heading 2
*Italics*
`code`
```
code block
```
