![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
please do not modify these first two lines of the .md file
use the same syntax to add pictures:
placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Calculating the burden of mutations using conservation scores

### Introduction

Another effect of being in a very small population size is accumulating deleterious mutations, or mutational load. With the increasing availability of whole genome sequences of various non-model species, we are now able to make a multi-species alignment of hundreds of animal species and look for regions that are conserved, or in another words, preferred by natural selection. When these regions experienced mutations, it means that the area that has been highly adapted and selected was disrupted, and this might not be good for the species' survival.

In this tutorial, we will calculate mutational load from a pre-determined site. What is needed in this tutorial is:
1. Determine the sites where we want to calculate mutational load.
2. Obtain the score reference file for those sites.
3. Getting the genotype likelihood of all genotypes from the pre-determined sites from all samples.
4. Getting the genotype probability of all genotypes (AA,AC,AT,AG,CC,CG,CT,GG,GT, TT) from all sites from the genotype likelihood.
5. Align the score reference to the genotype probability file.
6. Calculate the total mutational load within each sample.

Step 1 to 3 needs takes quite a while so it has been done for this tutorial. Step 3 especially takes ~1h per sample. The command to generate a genotype likelihood file is given below for your information.
```{bash glf, eval=FALSE}
FILE=input.bam
SITES=predeterminedSites.file # indexed in ANGSD

# Getting glf.gz file from a BAM
angsd -i $FILE -sites $SITES -out ${FILE%.bam} -minQ 20 -minMapQ 20 -remove_bads 1 -trim 5 -GL 2 -doGlf 4
```

Question: How do you decide which sites we use to calculate mutational load?

The conservation scores we will be using is the SIFT score. The SIFT score is a measure of deleteriousness of a mutation by how harmful is the change to the protein structure. Consequently, it is only available on coding sequence.

### Task 0: Prepping your working directory

> Exercise 0.1
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

Sometimes, we cannot be certain about the genotypes in a site because our sequencing coverage is low. In this situation, the preferred approach is to use the likelihood of a genotype which can be obtained with software tailored to analyse low coverage genome such as `angsd`. T

The genotype likelihood, however, only describes the likelihood of the genotypes relative to the best one. Have a look on one of the .glf file using `head` command. The resulting file should look like this:
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
This does not tell us how probable each genotypes are because the values are relative to the best one (0.0000000). We need a kind of "independent" value that tells us the genotype probability.

To do that, a Master student of ours, Deborah Greer, made a small custom script that can change genotype likelihood into genotype probability based on Bayesian Theorem together with Alberto Carmagnini.

You need to use the script to change the genotype likelihood file into genotype probability file using the provided script as below.
```{bash gpf, eval=FALSE}
python3 all_genotype_likelihoods_v3.py $FILE.glf $FILE.gpf
```

> Challenge: Make a one liner to do this simultaneously for all 18 .glf files!

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

> Food for thought: Do we have the same sites between .glf and .gpf? Are each of the samples having the same number of sites?

### Task 2: Aligning genotype probability file with conservation score reference

After obtaining the gpf file, we align the reference score that was premade (`SusScr11_107_sift_scores_chr1.bed`). Let's have a look on this file using `head()` command.
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
A bed file is a text file format used to store genomic regions as coordinates and associated annotations. The data are presented in the form of columns separated by spaces or tabs. Starting from the fourth column and subsequen columns, the content is usually variable, ranging from gene names, score, strand information, and many more (Read more (here)[https://en.wikipedia.org/wiki/BED_(file_format)]). In this tutorial, it is the ancestral allele of the site, and the deleteriousness score when the derived allele is A, C, G, and T (alphabetical order).

We would like to align this reference to our .gpf files for further processing downstream. We will use `bedtools intersect`.
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

Question: What is happening when we did the bedtools intersect? How many sites do we have now?

### Task 3: Calculating mutational load using custom script

The `_sift.bed` file is the input for another script that was initiated by Deborah Greer also, to calculate the total load per sample that are homozygous, heterozygous, or both. After the score has been aligned and made in the same file as the genotype probability score, we run the custom python script as follows:
```{bash load, eval=FALSE}
python mut_load_calculator_SIFT.py ${FILE%.gpf}_sift.bed ${FILE%.gpf}_0 0
```

Question (Advanced): Have a look on the mut_load_calculator_SIFT.py script. What do you think the argument '0' stands for?

The resulting files of this commands will be `*_sift_scores.txt` and `*_sift_het_scores.txt`. The first contain the homozygous load while the last contain the heterozygous load. We will discuss the content of these files when we are trying to plot the files.

### Task 4: Concatenating the results

As we have one file for each sample, we would like to have them in one file for easy plotting with R. 

What we will be doing is first concatenate the .txt file with same suffix in one file. For example, all files ending up with *_sift_scores.txt can be concatenated as follows.
```{bash cat, eval=FALSE}
cat $(ls result/RD*sift_scores.txt | sort -V) | sed 's/result\///' | sed 's/_chr1_sift.bed//' > babirusa_sift_scores.txt
```

> Exercise 4.1.
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

> Exercise 5.1.
>
> Merge the metadata with left_join() as in with previous session.

Then we plot the load as follows.
```{r plotLoad, eval=FALSE}
library(ggplot2)
ggplot(l)+
  geom_boxplot(aes(x=Region,y=homLoad))+
  labs(y="homozygous load")+
  theme_minimal()
```

Question: Which region has the highest homozygous load? Which has the lowest?

### Task 5 (Optional)

The other file, "*sift_het_scores.txt" contain more information such as the heterozygous ancestral and heterozygous derived, assuming that the ancestral allele we got are all homozygous.

When you have the time, repeat the commands from Task 3, adjusting the file names accordingly, and try to plot the result. The column headers are given below:
```{r hetFile, eval=FALSE}
c("SampleID","totalHet","totalHetAnc","totalHetDer",
                 "totalScoreHet","totalScoreHetAnc","totalScoreHetDer",
                 "hetLoad","hetAncLoad","hetDerLoad",
                 "totalLoad","totalLoad_hetAnc", "totalLoad_hetDer","hetFactor","scoreType")
```

Question: Is there any difference between heterozygous ancestral and derived? Would you trust the result?
