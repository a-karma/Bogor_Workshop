![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data procesing - Session 3
In this section we will highlight the necessary steps to identify genetic variations from whole-genome-sequencing data.
 
### Variant Calling Format (vcf)
VCFs are tab-separated text files representing the standard file format for storing genetic variation data. Each .vcf file should contain some header lines (starting with `#`) which provide metadata information describing the body of the file. The file body consists of 8 mandatory columns and an unlimited number of optional columns that may be used to record other information about the sample(s). When additional columns are used, the first optional column is used to describe the format of the data in the columns that follow (see table below)

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



##### Example of sub-section title 
- this is how you define bullet point
- in case you want to make more explicit 
- the series of steps required to accomplish a task 

See example below on how to format commands that the participants will have to run

```sh
conda activate Day_1
plink --bfile file_name --recode
```

you can instead use `this syntax` to highlight an in-line command, software name or something you think it's important

see below the syntax for tables:

| column A | column B |
| ------ | ------ |
| row 1a | row 1b |
| row 2a | row 2b |
| you can also leave cells blank | |

you can also use this env for exercises and tips:
> Exercise 1 
> 
> Modify the command above to ...

Or even:
> Best practice: never use spaces in file names

note the empty line in the first quoted env (it seems to make it look nicer on github) 
