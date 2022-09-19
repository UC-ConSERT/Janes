First README file


## 02 Aug 2022
Tara iti allignments failed: no space, but could be because running too many scripts at a time (have space)

# Tuturuatu
## Mosdepth/pcr dup removal
Using Markdup:cat
-
- Rerun mosdepth and python. Needs to be done within an environment. Error where mosdepth made file (global.dist.txt) that was not complete so was getting error      "ValueError: not enough values to unpack (expected 3, got 1)", as stopped part way through the making of the 3 columns in the global.dist.txt file. This is likely due to an error with running mosdepth outside an environment. Jana fixed by creating an environment with mosdepth within.
    - Ran 2_9 first for mosdepth and python. Need to run 3 for mosdepth and python? Or did I already do that? I think maybe at the same time as 2_9? I think I did. Evidence will be in if python works. How to open it?


If qualimap (.graphmap) doesn't look right for nodup files, this could be because I killed it part way through the first one (CR01).
In nodup stats, I think CR01.stats etc is just copies of dup.stats made earlier (woops). Accidentally deleted CR02.stats from this folder (but think this is a duplicated of dup.stats)

## Where I am up to
Need to make TLR bed file.
Finish and run counting snps script.
Import Molly's graphs scripts and run graphs. Find out how to run R online? Maybe use jupyter hub...or R online?
Change variant calling process to output vcf files? that can be gzipped later?
Filtering does weird things to file names by adding bcf.recode.bcf... I wonder if this is because it is expecting a vcf file (from Molly's previous steps in creating the concat file), which it then recodes to bcf..? Check manual, filter script includes "--recode bcf" or something...
# 19/9/22
Changed all of 5_... scripts from #!/bin/sh to #!/bin/bash -e . Haven't tested if they run yet.

Necessary programs:
- samtools
- bamtools
- bcftools
- etc...
- qualimap
- mosdepth

Information about data.
[Link text](link address)

# Heading 1
## Heading 2
*Italics*
