![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data procesing - Session 2
In this session you will learn how to align reads to a reference genome and how to manipulate bam files.

#### The SAM and BAM format
SAM stands for Sequence Alignment/Map format. It is a TAB-delimited text format consisting of a header
section, which is optional, and an alignment section. The BAM is the compressed binary version of a SAM file.
Header lines start with `@`, while alignment lines do not. Each alignment line has 11 mandatory fields for
essential alignment information such as mapping position, and variable number of optional fields for flexible
or aligner specific information. 

In your raw_data directory you should see a folder called CHR_10. 
This contains a few bam files that we have generated for you using the same procedure that you are learning today. 
Let's examine the first few line of the header:

```sh
samtools view -H RD10_chr10.bam | head -20
```


#### Indexing BAMS 
Indexing aims to achieve fast retrieval of alignments overlapping a specified region without going through
the whole alignments. BAM must be sorted by the reference ID and then the leftmost coordinate before
indexing.
