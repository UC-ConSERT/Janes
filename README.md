First README file


## 02 Aug 2022
Tara iti allignments failed: no space, but could be because running too many scripts at a time (have space)

Tuturuautu no dup doesn't exist. Markdup instead?
Using Markdup:
- Couldnt find nodup files after markdup? and before Marking duplicates for next file:
    [E::hts_open_format] Failed to open file "CR01_nodup.bam" : No such file or directory
    samtools stats: failed to open "CR01_nodup.bam": No such file or directory
- "Running calculating stats" failed. Index doesn't exist for the nodup files. Need to index using:
#Indexing the merged bam file
for file in ${sppdir}nodup_bam/*_nodup.bam
do
        base=$(basename $file _nodup.bam)
        echo "Indexing nodup bam file $base"
        samtools index -@ 16 -b ${sppdir}nodup_bam/*_nodup.bam
done
echo "Indexing merged bam files is complete"
- Rerun python at bottom wih correct file location

Using Markdup:
- Couldnt find nodup files after markdup? and before Marking duplicates for next file:
    [E::hts_open_format] Failed to open file "CR01_nodup.bam" : No such file or directory
    samtools stats: failed to open "CR01_nodup.bam": No such file or directory
- "Running calculating stats" failed. Index doesn't exist for the nodup files. Need to index using:
#Indexing the merged bam file
for file in ${sppdir}nodup_bam/*_nodup.bam
do
        base=$(basename $file _nodup.bam)
        echo "Indexing nodup bam file $base"
        samtools index -@ 16 -b ${sppdir}nodup_bam/*_nodup.bam
done
echo "Indexing merged bam files is complete"
- Rerun python at bottom wih correct file location



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
