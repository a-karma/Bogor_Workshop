3![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data procesing - Session 3
In this section we will highlight the necessary steps to identify genetic variations from whole-genome-sequencing data.

### Variant Calling procedure
Variant calling is a crucial step in analyzing next-generation sequencing (NGS) data, involving the identification of genetic variations between a reference genome and the genomes of individuals or populations. These variations, collectively known as variants, encompass single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variations. 

Several methods are employed for variant calling in NGS data analysis, each with its unique strengths and limitations. Here, we'll focus on alignment-based methods, which rely on aligning sequencing reads to a reference genome and identifying variants based on the discrepancies between the reads and the reference sequence. These methods include Samtools, BWA/GATK, and FreeBayes, with FreeBayes being the software we'll be using in this session.
 
### Variant Calling Format (vcf)
The VCF file format, also known as the Variant Call Format, is a standard text-based format for storing and exchanging information about genetic variations. It consists of tab-delimited text files with header and body sections.

The header section, identified by lines starting with #, provides metadata about the VCF file, such as the reference genome, variant caller,  and sample information. The body section contains detailed information about each variant in the sample(s). It consists of eight mandatory columns, which are as follows:

|Column Index| Name| Description |
|---|---|---|
|1|CHROM| The name of the sequence (typically a chromosome) on which the variation is being called.|
|2|POS| The 1-based position of the variation on the given sequence.|
|3|ID|The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons.|
|4|REF|The reference base (or bases in the case of an indel) at that position.|
|5|ALT|The list of alternative alleles at this position.|
|6|QUAL|A quality score associated with the inference of the given alleles.|
|7|FILTER|A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.|
|8|INFO|An extensible list of key-value pairs (fields) describing the variation.|
|9|FORMAT|An (optional) extensible list of fields for describing the samples.|
|+|SAMPLEs|For each (optional) sample described in the file, values are given for the fields listed in FORMAT| 

Now that we have familiarised ourselves with the file format let's have a look at how we can generate a vcf.

### Commands for calling variants with freeBayes
Today we going to use freeBayes, which is an easy to use, fast genetic variant detector designed to find SNPs (single-nucleotide polymorphisms) and small indels (insertions and deletions). You can read more about freeBayes here [github page](https://github.com/freebayes/freebayes)

The general usage of this variant calling software is as follow:
```sh
freebayes -f reference_genome -L list_of_bam_files -v output.vcf
```
Have a look at the various options and at the examples by running:
```sh
freebayes --help
```
Let's experiment with this software. First of all, we need to create a list of bam files the we want to include in our analysis.
```sh
ls ~/day2/raw_data/chr_10/*.bam > ~/day2/lists/bams_for_vcf.txt
```
Even with a small sample size like the one we have provided you with, the variant calling procedure is very computationally demanding. In order to complete the task in the time allowed for this session and to avoid server overload, we are going to examine only a small region of chromosome 10. The region can be specified with the `-r` flag as in the example below:
```sh
freebayes -f ~/day2/raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -r 10:0-1000000 -L ~/day2/lists/bams_for_vcf.txt > ~/day2/vcfs/babirusa.chr10.0to1mb.vcf
```
Here we have called variants only in a region of 1Mb at the beginning of chromosome 10. 

> `Exercise 1`
>
> run freeBayes on the same region but using a more stringent base quality filtering (-q 30) and mapping quality filtering (-m 20). 
> Then calculate how many variants have been included in the output and compare this value with the previous one. Did you expect this result?

Let's keep experimenting with the filtering options:

> `Exercise 2`
>
> run freeBayes on the same region but using a even more stringent filtering; -q 30 -m 20 -C 5.
>
> Have a look at the help message to understand the meaning of the `-C`
>
> Then calculate how many variants have been included in the output and compare this value with the previous ones. Did you expect this result?

Let's examine another region:

```sh
freebayes -f ~/day2/raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -r 10:1000000-2000000 -L ~/day2/lists/bams_for_vcf.txt > ~/day2/vcfs/babirusa.chr10.1to2mb.vcf
```
by changing the values after the `-r` flag we are moving along chromosome 10 and analysing the next region of 1Mb.

> `Question`: how many variants freeBayes have identified this time?

In principle we could write a script that splits the entire chromosome 10 into regions and pass this list of values to freebayes using a loop.
This goes beyond the scope of today but feel free to experiment with this if you have time.

We will now introduce a software that will become very useful for tomorrow's sessions that will allow us to concatenate vcf files: `vcftools`

First of all we need to compress and index the vcfs that we want to concatenate:
```sh
bgzip -c ~/day2/vcfs/babirusa.chr10.0to1mb.vcf > ~/day2/vcfs/babirusa.chr10.0to1mb.vcf.gz
bgzip -c ~/day2/vcfs/babirusa.chr10.1to2mb.vcf > ~/day2/vcfs/babirusa.chr10.1to2mb.vcf.gz
tabix -p vcf ~/day2/vcfs/babirusa.chr10.0to1mb.vcf.gz
tabix -p vcf ~/day2/vcfs/babirusa.chr10.1to2mb.vcf.gz
```
Then we can concatenate the two by running:
```sh
 vcf-concat ~/day2/vcfs/babirusa.chr10.0to1mb.vcf.gz ~/day2/vcfs/babirusa.chr10.1to2mb.vcf.gz | bgzip -c > ~/day2/vcfs/babirusa.chr10.0to2mb.vcf.gz
```

We will lear a lot more about vcfs and genetic diversity tomorrow.





