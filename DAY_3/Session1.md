![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation


## Day 3 - Population genetic analysis - Session 2
In this session you will learn how to measure genetic diversity and what can we learn from these metrics.
### Measuring Genetic Diversity 
We are going to explore the observables of genetic data a.k.a genetic diversity. We are going to measure three metrics of genetic diversity: Observed heterozygosiy (Ho), number of seggregating sites (S) and the average number of pairwise difference between sequences (Pi).
We will further try to explore if a population is evolving or is neutral using these metrics.

##### Heterozygosity 
- Heterozygosity is the probability of sampling two different alleles for a locus. Heterozygosity is an important predictor of survivability of populations. It helps rescue from consequences of recessive deleterious alleles for example.  
- We will use vcftools(https://vcftools.sourceforge.net) for this exercise 
- There are other methods as well, for example, RTG tools (https://github.com/RealTimeGenomics/rtg-tools) and ANGSD (http://www.popgen.dk/angsd/index.php/Heterozygosity)

First we will estimate the number of heterozygous genotypes at a locus. That is, the number of individuals that are heterozygous at a locus.
```sh
conda activate Day_3
# Calculate the number of heterozygous genotypes at a loci using vcftools
vcftools --vcf [input.vcf] --hardy --out [input]

```

Then we will estimate the proportion of individuals that are heterozygous at a locus. For this we need to calculate the ratio of number of individuals heterozygous at a locus and the total number of individuals genotyped at the locus.

Let's first visualize the output generated in the previous step. For this we will use the `less` command. This command allows us to read a file without printing it or opening the whole file. The `less -S` option will let us read the file without word wrapping. Use `ctrl+q` to exit the window.

```sh
# View the output  
less -S [input].hwe

```

Now, let's estimate the heterozygosity per locus. We will use the `tail -n +2` command. This will open the file without reading the first line which contains the headers.

Then, we print specific columns and perform the required mathematical operations for each column using the `awk` command. `awk` is a powerful command for performing column operations.

```sh
# Estimate heterozygosity per loci
tail -n +2 [input].hwe | awk -F "[\t/]" '{print $4, $3+$4+$5, $4/($3+$4+$5)}' | less

```

```sh
# Calculate the number of heterozygous genotypes for an individual using vcftools
vcftools --vcf [input.vcf] --het --out [input]

# View the output and calculate heterozygosity per individual
less -S [input].het
tail -n +2 [input].het | awk '{print $1, 1-($2/$4)}' | less

```
> What factors influence heterozysity?

> How do you estimate expected heterozygosity?

##### Number of seggregating sites (S)
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

##### Average number of pairwise difference between sequences (Pi)
- Number of difference between every pair of sequence in the dataset is estimated and averaged out.
- We will use vcftools (https://vcftools.sourceforge.net) for this exercise
- There are other methods as well, for example, ANGSD (http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests)


```sh

# Calculate Pi per loci using vcftools
vcftools --vcf [input.vcf] --site-pi --out [input]

# View the output and estimate average Pi
less -S [input].sites.pi

```
> What factors influence Pi?

> What is the relation between Pi and S?
