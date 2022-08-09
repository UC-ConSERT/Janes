First README file


## 02 Aug 2022
Tara iti allignments failed: no space, but could be because running too many scripts at a time (have space)

Tuturuautu no dup doesn't exist. Markdup instead?
Using Markdup:cat
- Couldnt find nodup files after markdup? and before Marking duplicates for next file:
    [E::hts_open_format] Failed to open file "CR01_nodup.bam" : No such file or directory
    samtools stats: failed to open "CR01_nodup.bam": No such file or directory
        Hopefully fixed by replacing base with bam in input.
- Rerun python at bottom wih correct file location

If qualimap (.graphmap) doesn't look right for nodup files, this could be because I killed it part way through the first one (CR01).
In nodup stats, I think CR01.stats etc is just copies of dup.stats made earlier (woops)


Necessary programs:
- samtools
- bamtools
- etc...
- qualimap
- mosdepth

Information about data.
[Link text](link address)

# Heading 1
## Heading 2
*Italics*
