![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
please do not modify these first two lines of the .md file
use the same syntax to add pictures:
placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Detecting Recent Inbreeding using Runs of Homozygosity

### Introduction
Runs of Homozygosity (ROHs) are stretches of DNA where both chromosomes are identical (homozygous). This happens after consanginous mating i.e. when an individual inherits the same allele from related inviduals. Long ROHs suggest recent inbreeding, potentially leading to 1) reduced genetic diversity, making populations less adaptable to change 2) increased risk of deleterious recessive traits to be expressed, impacting health and survival. 

ROHs are a powerful tool for population genomics, helping us understand the past, present, and future of species. This is why in the last decade, ROH analyses have become the state-of-the-art method for inbreeding assessment. There are many ways to detect segments of ROH in a whole genome sequence. We can do it straight from BAM file (e.g. using methods like ROHan: https://github.com/grenaud/ROHan), or from a VCF file (most common). The principle is to compute heterozgosity in "sliding windows" across the VCF file to assess for the preseence of long streches of homozygosity in the genome an individual.

Please speak to your instructor if you do not understand the concept of "sliding windows" along the genome.

In this tutorial, we will use the most commonly used method: `PLINK`. `PLINK` is probably the most widely use software for conducting SNP analyses, including computing ROHs. Have a look at the `PLINK` manual to get an idea of how comphrensive this tool is: https://www.cog-genomics.org/plink/  

In `PLINK`, the `--homozyg` function is used to perform ROH analyses and relies on several input settings.

### Task 0: Preparing your working directory

> Exercise 0.1
>
> Prepare your working directory for this session by following all the steps in Task 0 of Day 4 Session 1, but with data from /home/DATA/Day_4/Session_2

You should have a working directory in your home folder now named `day_4_inbreeding` containing a symbolic link of `babirusa_workshop_set.vcf.gz`, `babirusa_workshop_metadata.txt`, and 18 files of a ROH run results ending with `.hmmrohl.gz`. We will work with these results a bit later in the session.

### Task 1 Detecting ROH with PLINK

PLINK detect subsequent runs of homozygous genotypes in a genomic data sample using VCF data. To remind you, a VCF file contains variants from single or multiple sample.

> Exercise 1.1
>
> Have a look on the vcf file we are going to use and see how many samples we have. Use `bcftools view -h <filename> | tail -n 1`.

To detect homozygous segments in a VCF, we run a `plink` command with the `--homozyg` flag as follows:
```{bash plink, eval=FALSE}
plink --vcf babirusa_workshop_set.vcf.gz \
  --homozyg \
  --out babirusa_workshop_set_PLINK_A \
```

The output files of the command will be a set of text files with the prefix declared in the --out argument ending with suffix ".hom.*". The ".log" file will contain all the options and the command 
```sh
ls -lh babirusa_workshop_set_PLINK_A.*
 539K Dec  4 13:43 babirusa_workshop_set_PLINK_A.hom
  821 Dec  4 13:43 babirusa_workshop_set_PLINK_A.hom.indiv
 161M Dec  4 13:43 babirusa_workshop_set_PLINK_A.hom.summary
 1.2K Dec  4 13:43 babirusa_workshop_set_PLINK_A.log
  170 Dec  4 13:43 babirusa_workshop_set_PLINK_A.nosex
```

The detail of the ROH segments found can be found in the ".hom" file.
 ```sh
 head babirusa_workshop_set_PLINK_A.hom
 FID  IID      PHE  CHR SNP1 SNP2         POS1         POS2         KB     NSNP  DENSITY     PHOM     PHET
 RD1  RD1   -9.000    1    .    .      1198144      2953307   1755.164     1445    1.215    0.997    0.003
 RD1  RD1   -9.000    1    .    .      3585147      5219183   1634.037     2885    0.566    0.992    0.008
 RD1  RD1   -9.000    1    .    .      7043994      8521042   1477.049     2199    0.672    0.995    0.005
 RD1  RD1   -9.000    1    .    .     10434224     11627525   1193.302     2262    0.528    0.994    0.006
 RD1  RD1   -9.000    1    .    .     12676447     13694281   1017.835     2010    0.506    0.995    0.005
 RD1  RD1   -9.000    1    .    .     13697923     14982005   1284.083     2349    0.547    0.995    0.005
 RD1  RD1   -9.000    1    .    .     17237019     18584878   1347.860     2418    0.557    0.995    0.005
 RD1  RD1   -9.000    1    .    .     24939490     26049512   1110.023     2045    0.543    0.992    0.008
 RD1  RD1   -9.000    1    .    .     28662930     29897745   1234.816     2197    0.562    0.995    0.005
 ```
You can immediately recognise that this output file gives sample names, the start and end position of the ROH segment length detected, which chromosome it is, the number of SNP representing the segment. Read more about the output details in this [link](https://www.cog-genomics.org/plink/1.9/formats#hom).

> Exercise 1.2
>
> Have a look on the content of the "`.hom.indiv`" output. It will tell you the summary of ROH segment distribution per sample. How many samples in our sample set has ROH?

The power of ROH is that the segment length is inversely correlated with the time of inbreeding event. The longer the segment is, the more recent is the inbreeding event. However, we cannot see this easily from the output file. We need to download the output files and plot the results. Before we do that, we need to prepare a R working directory in our local computer.

### Task 2 Downloading and preparing PLINK results for plotting in R

We will use R in R Studio to visualise our ROH results. To make our working directory, open R Studio, choose File > New Project and in the `New Project Wizard` choose `Create New R Directory > New Project`. Type in your directory name, such as "Inbreeding_Analysis" and choose your preferred path within your local computer. If you are happy with the directory name and location, click `Create Project`. A new R Studio session will appear.

To start our data visualization project, choose `File > New File > R Script`. Save the R Script from the start to avoid forgetting to save your commands. For example, save it as "`01_plot_PLINK.R`" in our R project working directory.

Then, go to your command line interface and go to the directory of your newly made R project. Make a new directory called "`input`" and download your PLINK results there. Do not forget to also download the accompanying metadata you have sym-link-ed. For example:
```sh
mkdir input
cd input
sftp -i <path_to_identity_file> <username>@138.246.238.65
> get /home/<username>/day_4_inbreeding/babirusa_workshop_set_PLINK_A.hom* .
> get /home/<username>/day_4_inbreeding/babirusa_workshop_metadata.txt
```
Do not download the `.vcf` file as it will be veru heavy.

> Tips: You can check whether you are on the correct directory by looking at the content. If it contains a file with the directory name with .Rproj suffix, you are on the right place.

Afterwards, we go back to our R Studio and open our "`01_plot_PLINK.R`".

To plot the results, we need to first read the PLINK results and the accompanying metadata.
```{r readHom}
library(dplyr)
p<-read.table("input/babirusa_workshop_set_PLINK_A.hom", header=T)
m<-read.table("input/babirusa_workshop_metadata.txt",header=T)
```

To check whether our dataset were properly read into the assigned object `p` and `m`, we use the R command `head()`.
```{r}
head(p)
head(m)
```

We would like to plot while including the metadata. To do this, we will use `left_join()` function from the package `dplyr` to join the two.
```
r<-left_join(m,p,by=c("Sample"="IID"))
```
The argument `by=c("Sample"="IID")` allows you to merge similar columns. Look at the result using `head()` is the merging successful? Pay attention in the observation count of the objects.

Now we have a dataframe object `r` containing ROH segment length details.

### Task 3 Visualising ROH results using `ggplot2` and interpreting it

To plot the PLINK results, we will use an R package `ggplot2`. This package uses a different plotting grammar than base R plotting functionalities. In principle, you need to choose a type of plot, e.g. `geom_bar()` or `geom_point()`, then determine the aesthetics inside the function, i.e. the elements of the graph like x axis, y axis, how is it colored, etc, and run it with `ggplot()` function. read more about what `ggplot2` can do (here)[https://r-graph-gallery.com/ggplot2-package.html].

So, our first plot will be looking at how ROH segment distributed along chromosome 1. We will use `geom_segment()` here that could plot segment of ROH by specifying the start and end of the x and y coordinates.
```{r plotSeg}
library(ggplot2)
r %>%
  filter(CHR==1) %>%
  ggplot()+
  geom_segment(aes(x=POS1/1e6, xend=POS2/1e6,
                   y=Sample, yend=Sample, color=Region), size=3)+
  labs(x="Chromosome 1 (Mbp)",y="Sample ID")+
  theme_minimal()
```
Note the `r %>% filter(CHR==1) %>%` before the `ggplot()` function. As ggplot2 is part of the `tidyverse` package set, it can be run with the tidyverse piping grammar (`%>%`) to direct a function/object to another function, much like bash pipe (`|`). Using this functionality, we are asking R to only read chromosome 1 from the `r` object using `filter()` function before plotting the segment.

> Exercise 3.1.
>
> Plot the ROH segment in chromosome 2 by modifying the `filter()` command. Is the sample with the most abundant ROH still the same as in chromosome 1?

The segment plot is a nice way of looking whether ROH has been happening at the same exact location in all samples, which could be an indication of selection within the population. Despite the coloring per region, it is still quite difficult to discern which region has the most ROH.

To know which region has the most ROH, we need to summarise the ROH length per sample. This can be done with tidyverse functionality as follows.
```{r group}
rSROH<-r %>% group_by(Sample) %>% summarise(SROH=sum(KB))
```
Open rSROH and see what is inside. This is the summary of the total ROH segment per sample (SROH) we obtained from PLINK.

> Exercise 3.2
>
> We have this summary on the .hom.indiv. Can you check whether we have the correct sum per sample?
> Repeat the above command to get the mean length and also check.

> Exercise 3.3.
>
> Merge `rSROH` with `r` using `left_join()` into a new object calles `r2`. Have a look first on the name of the columns that are the same in each object, i.e. defining the SampleID, to use the `by` argument (See Task 2).

<details open>
<summary>Exercise 3.3. Solution</summary>
<br>
```{r}
r2<-left_join(r,rSROH,by=c("Sample"="Sample"))
```
</details>

To see which population has the most ROH, we will group the SROH values based on the region and made a boxplot.
```{r}
ggplot(r2)+
  geom_boxplot(aes(x=Region,y=SROH/1e6,color=Region))+
  labs(x="region",y="Total ROH amount (Mbp)")+
  theme_minimal()
```
Which population has the highest amount of inbreeding?

> Exercise 3.4.
>
> The inbreeding coefficient (FROH) is the total sum of ROH segment divided by the entire genome length in which we run the inbreeding analysis. As we map to the domestic pig genome and ran this on autosomal genome only, the total here is 2265774640 base pairs. Calculate the FROH and plot the result. Is it still looking the same with the SROH plot?

This plot, however, does not tell us how many are the short segments and long segments so that we can see how recent the inbreeding has been. To do that, we will use `geom_freqpoly()` to show the distribution of ROH segment lengths in a sample.
```{r}
ggplot(r2)+
  geom_freqpoly(aes(x=KB,color=Region, group=Sample),
                size=2, alpha=0.5)+
  xlim(1000,10000)+
  theme_minimal()
```
Now we know that most segments are short, and there is no distinct difference in segment length distribution between the different regions. This means that the extent of inbreeding is almost natural, where it is mostly broken down by recombination in the following generations and resulting in short segments.

Let's go back to the server and do Exercise 3.5.

> Exercise 3.5.
>
> We have run the default PLINK option of the `--homozyg` argument. The entire default parameter setting was based on human genome. A babirusa genome, however, was different from human genome. It may have different mutation rate, recombination rate, linked genes, etc.  that was based on human genome. However, we did not know any of this for the babirusa genome. We will approach it using the pig genomes that has been run by previous studies (reviewed by Meyermans et al. 2020), i.e., we will run a second PLINK with output suffix `_PLINK_B` with minimum SNP count 20 (`--homozyg-snp 20`), minimal ROH segment length 10 kb (`--homozyg-kb 10`), SNP density 1 per 1 Mbp (`--homozyg-density 1000`), and scanning window consists of 20 SNP (`--homozyg-window-snp 20`) with 0.25 portion of the overlapping windows must be called homozygous to define any given SNP as 'in a homozygous segment' (`--homozyg-window-thershold 0.25`).
>
> Read the documentation for PLINK `--homozyg` (here)[https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo] for more information about the parameters.
>
> Compare the result of Exercise 2.4 with the result from the default. Which has more ROHs per sample? Which set of parameters do you trust more?

<!---
I just added this, might not be a good idea? Just intense plotting session.
--->

### Task 4: Comparing observational method with model-based method

The many input parameters required by observational method such as PLINK makes ROH results volatile. Another way to detect ROH is to use a model-based method, i.e., given a certain model on how to detect a homozygous genotypes, how likely is a segment a ROH. Such model-based method can be really long to run (~6 hours with 18 threads per sample!), so we have run the model for you and we will focus on analysing the results instead. 

The model was run using ROHan (Renaud et al., 2018) with only the mutation rate in ROH parameter set to 0.0001. We will discuss more about this method during the tutorial's wrap up session.

If you have done Exercise 0.1, these results should have been in your `day_4_inbreeding` directory. Let's have a look on one of the results using `zcat` and `head`.
```sh
#ROH_ID	CHROM	BEGIN	END	ROH_LENGTH	VALIDATED_SITES
1	1	67000001	68000000	1000000	726880
2	1	94000001	96000000	2000000	1420468
3	1	117000001	118000000	1000000	688128
4	1	122000001	123000000	1000000	637574
5	1	126000001	128000000	2000000	1238786
6	1	173000001	181000000	8000000	5120675
7	1	192000001	195000000	3000000	2066042
8	1	208000001	209000000	1000000	767864
9	1	214000001	217000000	3000000	1924849
```
This looks similar to our PLINK `.hom` results. This means we can plot them similarly.

Download these results into your R project so that it is on the same directory with your PLINK results. Then, read the files into R using `read.tables()`, add a column containing Sample ID, and concatenate all files into one file to make a similar dataframe as we have with PLINK where it contain all details of the ROH of all samples.

For example:
```{r}
s1<-read.table("input/RD1_mdup_1e-4.mid.hmmrohl.gz")
colnames(s1)<-c("rohID","chr","start","end","length","validatedSites")
s1$SampleID<-"RD1"
```

> Challenge: How can you make your task faster with R? Hint: make a list of the file paths.

If you do not feel challenged, feel free to repeat the commands manually or run the following code (Warning: it needs you to install an optional package `janitor` and `gtools`).
```r
library(gtools)
library(janitor)

s<-list.files(path="input/",pattern="*.mid.hmmrohl.gz",full.names = T)

rohan_df<-function(temp) {
  temp<-mixedsort(temp) 
  list<-lapply(temp,read.delim)
  ID_list<-unlist(regmatches(temp,gregexpr("[A-Z]+[0-9]+",temp)))
  for(i in 1:length(list)){
    list[[i]]$SampleID <- rep(ID_list[i],nrow(list[[i]]))
  }
  df<-do.call(rbind,list)
  return(df)
}

s_df<-rohan_df(s)
s_df<-clean_names(s_df, "lower_camel")
```

The resulting file should look like more or less like this:
```
> head(s_df)
  xRohId chrom    begin      end rohLength validatedSites sampleId
1      1     1 39000001 46000000   7000000        4918042      RD1
2      2     1 57000001 58000000   1000000         644816      RD1
3      3     1 59000001 65000000   6000000        4217317      RD1
4      4     1 68000001 69000000   1000000         752473      RD1
5      5     1 71000001 78000000   7000000        4562652      RD1
6      6     1 87000001 92000000   5000000        3000020      RD1
```

> Exercise 4.1.
>
> Plot the segment distribution, the FROH boxplot, and the frequency of segment length classes as you have done with PLINK in Task 3. How similar it is with the PLINK results?

Depending on the parameter and tools you used, you can have different ROH distribution inferred per population. This is why understanding the genomic architecture of the species in full is important for accurate inference of inbreeding.

### Challenge (Optional): Calculating the number generations since the last inbreeding

The useful part of directly observing recent inbreeding is that you will be able to know when is the recent inbreeding given the rate of recombination rate within the species' genome. As the recombination rate of the babirusa genome is unknown, we will work with 1 cM ~ 1 Mb.

Here, we will make an additional column containing the number of generations represented by each ROH segment.
```
library(dplyr)
s_df<-s_df %>% mutate(generations=100/(rohLength*2))
```
Note that the formula in the new column `generations` is 100 divided by rohLength multiplied twice. This comes from the expectation of ROH segment length to follow an exponential distribution with mean equals to 100 / ( 2 * g * c ), with g is generation and c is recombination rate.

We can plot the `generation` using `geom_histogram()` and see how the number of generations since inbreeding differs between different population.

> Challenge: How will the distribution of generation time change if the recombination rate in pig genome is known to be 0.76 cM/Mb ? Modify the above command accordingly and plot the results using `geom_histogram()`.
