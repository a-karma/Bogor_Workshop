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
This contains a five bam files that we have generated for you using the same procedure that you are learning today. 
Let’s have a look at a BAM file and play around with a versatile tool to manipulate this file
format: ![samtools](http://www.htslib.org/).

First try to type:
```sh
less ~/day2/raw_data/CHR_10/RD10_chr10.bam
and as usual type Q to close it
```

This is how a binary file looks like, not very user friendly or human readable right?

Let’s find a better way to visualize it. We should start by examining the header section:

```sh
samtools view -H ~/day2/raw_data/CHR_10/RD10_chr10.bam | head -20
```
The `view` command of samtools allows you to read, print, and convert SaM/BAM/CRAM files. The `-H` flags outputs only the header.
As you can see, there are 19 lines starting with `@SQ`. Each of these lines corresponds to a chromosome of the reference genome and the value you read after `LN:` corresponds to the chromosome length in base pairs (bp).  

Now let’s have a look at an alignment line. We can use again the view command:
```sh
samtools view ~/day2/raw_data/CHR_10/RD10_chr10.bam | head -1
```
As you can see there are two fields occupying the largest portion of the line: one is the sequence (10 th field) and the other one is the Phred-scaled base quality (11 th field).
Those correspond respectively to the second and 4th line of the fastq file used to generate this bam. This is definitely better than using less but what if we wanted to visualize a specific region of this chromosome?

You may have noticed that in the CHR_10 folder there are also 5 files with the .bai extension. 
Those are the indexes of the bams and they allow fast access to specific regions without going through the whole alignment.
thanks to the indexing procedure we can use another samtools command and look at any given portion of chromosome 10:
```sh
samtools tview ~/day2/raw_data/CHR_10/RD10_chr10.bam -p 10:21100
```

We are looking at the alignment on chromosome 10 starting at the position 21100 from
the beginning of the chromosome (the `-p` flag stands for position). The second line is full
of Ns because we haven’t provide any reference genome in a fasta format to this samtools commad. 
The third line is the consensus sequence while each of the following line represent a specific read that was mapped to that region of the reference genome. 
As you can see, reads are displayed in different colours which represent quality scores. Just type `?` to look at the visualization options. Press for example `n` or `c` and see what happens.

Now that we have familiarised ourselves with this file format, let's try to generate some new bams!

#### Aligning reads
Aligning short-read sequences is the foundational step to most genomic and transcriptomic analyses. 
There are several tools to perform this task such as Bowtie2, BWA, HISAT2, MUMmer4, STAR, and TopHat2 just to name a few.
Here we are going to focus on the Burrows Wheeler Aligner (BWA) but the procedure can be generalised to any aligner with minor modifications.



#### Indexing BAMS 
Indexing aims to achieve fast retrieval of alignments overlapping a specified region without going through
the whole alignments. BAM must be sorted by the reference ID and then the leftmost coordinate before
indexing.
