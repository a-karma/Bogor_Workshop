![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
please do not modify these first two lines of the .md file
use the same syntax to add pictures:
placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Calculating the burden of mutations using conservation scores

### Introduction

In small populations, the accumulation of harmful mutations accelerates due to the combined effects of genetic drift and inbreeding. This accumulation, known as mutational load, compromises overall health, adaptability, and reproductive success, ultimately increasing the extinction risk for these vulnerable populations. This effect is further exacerbated by inbreeding-induced reductions in heterozygosity, which can lead to the expression of deleterious recessive traits.

Genomic data provides a powerful tool for assessing mutational load at the individual level. By analyzing multi-species alignments encompassing hundreds of animal species, we can identify conserved regions, representing areas of the genome that have been maintained in the same way (i.e. did not accumulate mutations) for millions of years due to natural selection. Mutations within these conserved regions potentially disrupt highly adapted and selected sequences and are therefore likely to be deleterious.

In this tutorial, we will calculate mutational load from a list of pre-determined sites. The general pipeline runs from BAM file as follows:
1. Determine the sites which we want to include in our mutational load analysis.
2. Obtain the deleterious score these sites from a database.
3. Obtain the genotype likelihood for all 10 possible genotypes at the same sites for all samples.
4. Intersect files in step 1 and 2; i.e. determine whether the genotype possessed by an individual at a specific position in the genome is likely to be deleterious or not.
5. Use the result of 4 to calculate the total, genome-wide, mutational load for each sample.

Step 1 to 3 are very slow so we will not go through these steps during this tutorial. In fact, step 3 takes ~1h per sample. The command, in `angsd`, which we used to generate a genotype likelihood file is given below for your information:
```{bash glf, eval=FALSE}
FILE=input.bam
SITES=predeterminedSites.file # indexed in ANGSD

# Getting glf.gz file from a BAM
angsd -i $FILE -sites $SITES -out ${FILE%.bam} -minQ 20 -minMapQ 20 -remove_bads 1 -trim 5 -GL 2 -doGlf 4
```

The conservation scores we will be using is the SIFT score. The SIFT score is a measure of amino acid similarity, which indicates the potential impact of a mutation on protein function. Ranging from 0.0 to 1.0, the score signifies the likelihood of a mutation affecting protein function. Scores closer to 0.0 suggest a higher probability of deleterious effects, while scores closer to 1.0 indicate a lower probability of functional disruption. Variants with SIFT scores between 0.0 and 0.05 are generally considered deleterious. SIFT scores can only be computed for mutation in coding region of the genome. 

> Food for thought: How do you decide which sites we use to calculate mutational load?

Conservation scores such as SIFT is only one among the many ways of characterizing deleterious variations empirically. Other ways are looking directly at amino acids or RNA transcripts. Read more about detecting deleterious variations in natural populations [here](https://doi.org/10.1146/annurev-animal-080522-093311). 

### Task 0: Prepping your working directory

> Exercise 1
>
> Prepare your working directory for this session by following all the steps in Task 0 of Day 4 Session 1, but with data from /home/DATA/Day_4/Session_3.

By the end of this exercise, you should have separate working directory for this session containing 18 files with .glf suffix, two python scripts, one bed file containing SIFT reference scores and one metadata file.
```
~/day_4_mutationLoad$ ls
RD10_chr1.glf  RD1_chr1.glf   RD3_chr1.glf   RD56_chr1.glf  RD61_chr1.glf  RD7_chr1.glf                       all_genotype_likelihoods_v3.py
RD16_chr1.glf  RD20_chr1.glf  RD44_chr1.glf  RD59_chr1.glf  RD64_chr1.glf  RD8_chr1.glf                       babirusa_workshop_metadata.txt
RD17_chr1.glf  RD2_chr1.glf   RD53_chr1.glf  RD60_chr1.glf  RD71_chr1.glf  SusScr11_107_sift_scores_chr1.bed  mut_load_calculator_SIFT.py
```

### Task 1: Obtaining genotype probability from genotype likelihood

When calculating mutational load we want to make sure we got the genotype right. Getting the genotype wrong at a specific position could mean we are expecting an individual to posses a highly harmful mutation while it does not, which can badly affect our calculations. Have a read on the impact of low depth to genotype calling [here](https://www.researchgate.net/publication/353052622_A_beginner's_guide_to_low-coverage_whole_genome_sequencing_for_population_genomics) or [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7115901/). Sabhrina, Alberto, and Deborah have designed a method to do this, which is what we are going to teach you today. There are other methods out there too, for more information please have a look at this review when you are done: 

To take into account genotyping uncertainity the method uses angsd to calculate genotype likelhood which outputs a file that looks like this: 
```{bash glf_view, eval=FALSE}
1   14852   -57.707866    -57.707866    -57.707866    -4.158484   -57.707866	-57.707866    -4.158484   -57.707866    -4.158484   0.000000
1   15068   -67.325843    -4.851565   -67.325843    -67.325843    0.000000    -4.851565   -4.851565   -67.325843    -67.325843    -67.325843
1   15164   -67.325843    -4.851565   -67.325843    -67.325843    0.000000    -4.851565   -4.851565   -67.325843    -67.325843    -67.325843
1   23911   -57.707866    -57.707866    -4.158484   -57.707866    -57.707866    -4.158484   -57.707866    0.000000    -4.158484   -57.707866
1   38704   -48.089888    -48.089888    -48.089888    -3.465403   -48.089888    -48.089888    -3.465403   -48.089888    -3.465403   0.000000
1   38730   -54.941796    -54.941796    -54.941796    -4.157494   -54.941796    -54.941796    -4.157494   -54.941796    -4.157494   0.000000
1   128269    -57.707866    -57.707866    -4.158484   -57.707866    -57.707866    -4.158484   -57.707866    0.000000    -4.158484   -57.707866
1   131076    0.000000    -2.079242   -2.079242   -2.079242   -28.853933    -28.853933    -28.853933    -28.853933    -28.853933    -28.853933
1   132504    -28.853933    -2.079242   -28.853933    -28.853933    0.000000    -2.079242   -2.079242   -28.853933    -28.853933    -28.853933
1   132540    -48.089888    -3.465403   -48.089888    -48.089888    0.000000    -3.465403   -3.465403   -48.089888    -48.089888    -48.089888
```
This is a tab separated file which give us the log likelihood of each 10 possible genotype at the site for one individual genome. We can use this value to calculate the probability that this individual possess any specific genotype using this data using Bayesian statistics. We have provided you a small python script that can do this. 

To run this script type:

```{bash gpf, eval=FALSE}
python3 all_genotype_likelihoods_v3.py $FILE.glf $FILE.gpf
```

> Challenge: write a one line script to automate this for all 18 .glf files!

Have a look at the resulting .gpf file using the `head` command. The resulting file should look like this:
```{bash gpf_view, eval=FALSE}
1   14851   14852   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0
1   15067   15068   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0
1   15163   15164   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0
1   23910   23911   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0
1   38703   38704   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0
1   38729   38730   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0
1   128268    128269    0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0
1   131075    131076    0.99998   1e-05   1e-05   1e-05   0.0   0.0   0.0   0.0    0.0   0.0
1   132503    132504    0.0   1e-05   0.0   0.0   0.99998   1e-05   1e-05   0.0    0.0   0.0
1   132539    132540    0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0
```

> Food for thought: can you spot the formatting difference between the .glf and .gpf files?

### Task 2: Aligning genotype probability file with conservation score reference

We now need to intersect SIFT scores and the gpf file. The SIFT score based on the pig genome are found in this file `SusScr11_107_sift_scores_chr1.bed`. Let's have a look on this file using `head()` command.
```sh
head SusScr11_107_sift_scores_chr1.bed 
1	204560	204561	A	1	0	0	0
1	204561	204562	C	0	1	0	1
1	204562	204563	A	1	1	0	1
1	204563	204564	C	0	1	0	1
1	204565	204566	A	1	0	0	0
1	204566	204567	C	1	1	0	0
1	204567	204568	C	0	1	0	1
1	204568	204569	A	1	0	0	1
1	204569	204570	C	1	1	0	1
1	204571	204572	T	0	0	0	1
```
This is again a `BED` file as in day one. You can see how we have customise this bed file - we kept its basic structure: tab separated 0-based coordinates but from field 4 we have added our own informations: the ancestral (non-deleterious) state of the allele, and the SIFT score of all possible homozygous genotpye, A, C, G, T at this position. 

We would like to intersect this file with our .gpf files to merge both genotype probabilities and their deleterious (SIFT) score in the same file. To do this we will use `bedtools intersect` (on day 1)
```{bash bedtools, eval=FALSE}
bedtools intersect -b $FILE -a SusScr11_107_sift_scores_chr1.bed -wb | cut -f 1-8,12-21 > ${FILE%.gpf}_sift.bed
```

The resulting *_sift.bed file should look as follows:
```{bash bed_view, eval=FALSE}
1	513572	513573	T   1   0.02    1   1   0.0   0.0   1e-05   0.0   0.0   1e-05   0.0   0.99998   1e-05   0.0
1	513593	513594	T   0.02    1   1   1   0.0   0.0   0.0   4e-05   0.0   0.0    4e-05    0.0   4e-05   0.99988
1	515134	515135	T   1   1   0.03    1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0
1	516060	516061	A   1   1   0.01    1   0.0   0.0   0.0   0.0   1.0   0.0   0.0   0.0   0.0   0.0
1	525966	525967	T   0   1   0   1   0.0   0.0   0.0   0.0   1.0   0.0   0.0    0.0    0.0   0.0
1	531650	531651	T   0.03    1   1   1   0.0   1e-05   0.0   0.0   0.99998    1e-05   1e-05   0.0   0.0   0.0
1	532250	532251	A   1   1   0.03    1   0.0   0.0   1e-05   0.0   0.0   1e-05	0.0   0.99998   1e-05   0.0
1	535817	535818	T   0   1   0.01    1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0
1	541103	541104	C   1   1   0   0   1.0   0.0   0.0   0.0   0.0   0.0   0.0    0.0    0.0   0.0
1	546672	546673	A   1   0.02    1   1   0.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0   0.0   0.0
```
Note that we have merged the content of the two files side by side here.

Question: Why did we use a cut command after the intersect?
>Hint: run the same command as above without the `cut` part instead use `| less -S`

### Task 3: Calculating mutational load using custom script

The `_sift.bed` file is the input for another script that we will use to calculate the total load per sample. The python script as follows:
```{bash load, eval=FALSE}
python mut_load_calculator_SIFT.py ${FILE%.gpf}_sift.bed ${FILE%.gpf}_0 0
```

Question (Advanced): Have a look on the command we use to run `mut_load_calculator_SIFT.py` script above. What do you think the argument '0' means? Hint: have a look inside the python script.

The script will output two files 1) `*_sift_scores.txt` and 2) `*_sift_het_scores.txt`. The first contain the homozygous load while the last contain the heterozygous load. We will discuss the content of these files later in the session when are plotting the results.

### Task 4: Concatenating the results

As we have one file for each sample, we would like to have them in one file to make it easier to plot the results with R. 

What we will be doing is first concatenate the .txt file with same suffix in one file. For example, all files ending up with *_sift_scores.txt can be concatenated as follows.
```{bash cat, eval=FALSE}
cat $(ls result/RD*sift_scores.txt | sort -V) | sed 's/result\///' | sed 's/_chr1_sift.bed//' > babirusa_sift_scores.txt
```

> Question: What do you think the `sed` command does? Have a look on the resulting file when you remove one of the sed command and compare it with a file which has the entire commands intact.

> Exercise 2
>
> Create a new R project for plotting the result of the mutational load run and download the results into a 'input' directory using sftp.

### Task 5: Visualizing homozygous load per population

We will read the resulting file as follows.
```{r readFile, eval=FALSE}
s<-read.table("input/babirusa_sift_scores.txt")
colnames(s)<-c("SampleID","totalPositions",
                 "totalHomozygousProbability",
                 "totalAncestralProbability",
                 "totalScore","totalScoreHet","totalScoreHom",
                 "transversionLoad","hetLoad","homLoad",
                 "totalLoad","hetFactor","scoreType")
head(s)
```
> Question: what do you think the colnames() function in R is for? Why do we need it?

> Exercise 3
>
> Merge the metadata with left_join() as in with previous session.

Then we plot a boxplot:
```{r plotLoad, eval=FALSE}
library(ggplot2)
ggplot(l)+
  geom_boxplot(aes(x=Region,y=homLoad))+
  labs(y="homozygous load")+
  theme_minimal()
```

Question: Which population of babirusa has the highest homozygous load? Which has the lowest?

Based on what you see with the PSMC and ROH distribution of the babirusa population from four different areas, would you think the abundance of homozygous load make sense here? 

### Challenge (Optional)

The other file, "`*sift_het_scores.txt`" contain more information on the heterozygous load. This file is automatically generated as part of the script `mut_load_calculator_SIFT.py` when we ran it on Task 3.

When you have the time, repeat the commands from Task 4 with the "`*sift_het_scores.txt`" files until the plotting. You will only need the `hetLoad` part to plot the heterozygous load.

The column headers are given below:
```{r hetFile, eval=FALSE}
c("SampleID","totalHet","totalHetAnc","totalHetDer",
                 "totalScoreHet","totalScoreHetAnc","totalScoreHetDer",
                 "hetLoad","hetAncLoad","hetDerLoad",
                 "totalLoad","totalLoad_hetAnc", "totalLoad_hetDer","hetFactor","scoreType")
```

> Question: Is there any difference between heterozygous ancestral and derived? Would you trust the result?
