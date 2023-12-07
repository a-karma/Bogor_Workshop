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


You can read more about freeBayes on the project [github page](https://github.com/freebayes/freebayes)

the general usage of this variant calling software is as follow:
```sh
freebayes -f reference_genome -L list_of_bam_files > output.vcf
```
Have a look at the varius options and at the examples by running:
```sh
freebayes --help
```

```sh
freebayes -f ~/day2/raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -q 20 -m 10 -r 1 -L bam_list.txt > babirusa.chr10.r1.vcf
```

> `Exercise 1`
>
> run freeBayes on the same region but using a more stringent quality filtering (-q 30) and a
> then calculate how many variants have been included in the output 
