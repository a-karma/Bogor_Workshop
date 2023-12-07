![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data procesing - Session 3
In this section we will highlight the necessary steps to identify genetic variations from whole-genome-sequencing data.

### Variant Calling procedure
Variant calling is the process of identifying differences between a reference genome and the genome of an individual or a population of individuals. These differences, known as variants, include single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variations. Variant calling is an important step in NGS data analysis, as it provides insights into the genetic variation that underlies phenotypic differences between individuals and populations.

There are several methods used for variant calling in NGS data analysis, each with its advantages and disadvantages. Here we will focus on alignment-based methods which involve mapping the sequencing reads to a reference genome and identifying variants based on the differences between the reads and the reference genome. These methods include Samtools, BWA/GATK, and FreeBayes (which is the softare that we are going to use in this session)
 
### Variant Calling Format (vcf)
The standard file format for storing information about genetic variation is termed VCF. These are tab-separated text files contain some header lines (starting with `#`) which provide metadata information followed by the body of the file. The file body consists of 8 mandatory columns and an unlimited number of optional columns that may be used to record other information about the sample(s). When additional columns are used, the first optional column is used to describe the format of the data in the columns that follow (see table below)

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

Now that we have familiarised ourselves with the file format let's have a look at the algorithmic procedure to generate these files.

### Commands for calling variants with freeBayes
freeBayes is a fast genetic variant detector designed to find SNPs (single-nucleotide polymorphisms) and small indels (insertions and deletions). You can read more about freeBayes on the project [github page](https://github.com/freebayes/freebayes)

The general usage of this variant calling software is as follow:
```sh
freebayes -f reference_genome -L list_of_bam_files -v output.vcf
```
Have a look at the varius options and at the examples by running:
```sh
freebayes --help
```
Let's experiment with this software. First of all, we need to create a list of bam files the we want to include in our analysis.
```sh
ls ~/day2/raw_data/CHR_10/*.bam > ~/day2/lists/bams_for_vcf.txt
```
Even with a small sample size like the one we have provided you with, the variant calling procedure is very computationally demanding. In order to complete the task in the time allowed for this session and to avoid server overload, we are going to examine only a small region of chromosome 10. The region can be specified with the `-r` flag as in the example below:
```sh
freebayes -f ~/day2/raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -r 10:0-1000000 -L  ~/day2/lists/bams_for_vcf.txt >  ~/day2/vcfs/babirusa.chr10.0to1mb.vcf
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
freebayes -f ~/day2/raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -r 10:1000000-2000000 -L  ~/day2/lists/bams_for_vcf.txt >  ~/day2/vcfs/babirusa.chr10.1to2mb.vcf
```
by changing the values after the `-r` flag we are moving along chromosome 10 and analysing the next region of 1Mb.

> `Question`: how many variants freeBayes have identified this time?

In principle we culd write a script that splits the entire chromosome 10 into regions and pass this list of values to freebayes using a loop.
This goes beyond the scope of today but feel free to experiment with this if you have time.

We will now introduce a software that will become very useful for tomorrow's sessions that will allow us to concatenate vcf files: `vcftools`



