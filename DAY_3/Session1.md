![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
please do not modify these first two lines of the .md file

use the same syntax to add pictures:

placeholder name within square brackets and ../IM/file_name.png within parentheses

## Measuring genetic diversity
In this session we are going to explore the observables of genetic data a.k.a genetic diversity. We are going to measure three metrics of genetic diversity: Observed heterozygosiy (Ho), number of seggregating sites (S) and the average number of pairwise difference between sequences (Pi).
We will further try to explore if a population is evolving or is neutral using these metrics.

### Measuring Genetic Diversity 
##### Heterozygosity 
- Heterozygosity is the probability of sampling two different alleles for a locus.
- We will use vcftools (https://vcftools.sourceforge.net) for this exercise 
- There are other methods as well, for example, RTG tools (https://github.com/RealTimeGenomics/rtg-tools) and ANGSD (http://www.popgen.dk/angsd/index.php/Heterozygosity)


```sh
conda activate Day_3
# Calculate the number of heterozygous genotypes at a loci using vcftools
vcftools --vcf [input.vcf] --hardy --out [input]

# View the output and estimate heterozygosity per loci
less -S [input].hwe
tail -n +2 [input].hwe | awk -F "[\t/]" '{print $4, $3+$4+$5, $4/($3+$4+$5)}' | less

# Calculate the number of heterozygous genotypes for an individual using vcftools
vcftools --vcf [input.vcf] --het --out [input]

# View the output and calculate heterozygosity per individual
less -S [input].het
tail -n +2 [input].het | awk '{print $1, 1-($2/$4)}' | less

```
> What factors influence heterozysity?

> How do you estimate expected heterozygosity?

##### Number of seggregating sites 
- Number of seggregating sites are the total number of polymorphic loci.
- We will use vcftools (https://vcftools.sourceforge.net) for this exercise 


```sh

# Calculate the number of polymorphic loci using vcftools
vcftools --vcf [input.vcf] --counts --out [input]

# View the output and estimate the number of seggreating sites
less -S [input].frq.count
tail -n +2 [input].frq.count | awk 'NF==6 {print $0}' | awk -F "[\t:]" '$6>0 && $8>0 {print $0}' | wc -l

```
> What factors influence number of seggregating sites?



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
