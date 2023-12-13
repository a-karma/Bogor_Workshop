![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
  
  <!---
  please do not modify these first two lines of the .md file
  use the same syntax to add pictures:
  placeholder name within square brackets and ../IM/file_name.png within parentheses
  --->
  
## Day 4 - Detecting Recent Inbreeding using Runs of Homozygosity - Session 2
  
### Introduction
Runs of Homozygosity (ROHs) are stretches of DNA where both chromosomes are identical (homozygous). This happens after consanginous mating i.e. when an individual inherits the same allele from related inviduals. Long ROHs suggest recent inbreeding, potentially leading to 1) reduced genetic diversity, making populations less adaptable to change 2) increased risk of deleterious recessive traits to be expressed, impacting health and survival. 

ROHs are a powerful tool for population genomics, helping us understand the past, present, and future of species. This is why in the last decade, ROH analyses have become the state-of-the-art method for inbreeding assessment. There are many ways to detect segments of ROH in a whole genome sequence. We can do it straight from BAM file (e.g. using methods like ROHan: https://github.com/grenaud/ROHan), or from a VCF file (most common). The principle is to compute heterozgosity in "sliding windows" across the VCF file to assess for the preseence of long streches of homozygosity in the genome an individual.
  
Please speak to your instructor if you do not understand the concept of "sliding windows" along the genome.
  
In this tutorial, we will use the most commonly used method: `PLINK`. `PLINK` is probably the most widely use software for conducting SNP analyses, including computing ROHs. Have a look at the `PLINK` manual to get an idea of how comphrensive this tool is: https://www.cog-genomics.org/plink/  
  
In `PLINK`, the `--homozyg` function is used to perform ROH analyses and relies on several input settings.
  
### Task 0: Preparing your working directory
  
> `Exercise 1`
>
> Prepare your working directory for this session by following all the steps in Task 0 of Day 4 Session 1.
>
> This time your project directory should be called `day_4_inbreeding` and the path to the input data for the link command should be `/home/DATA/Day_4/Session_2` (`ln -s /home/DATA/Day_4/Session_2/* .`).
  
You should have a working directory in your `/de/workshop_participants/your_user_ID/` folder now named `day_4_inbreeding` containing a symbolic link of `babirusa_workshop_set.vcf.gz`, `babirusa_workshop_metadata.txt`, and 18 files of a ROH run results ending with `.hmmrohl.gz`. We will work with these ROH run results a bit later in the session.
  
### Task 1: Detecting ROH with PLINK
  
First, we will detect runs of homozygous segments in a VCF using a software called `plink`. We run a plink command as follows:
```{bash plink, eval=FALSE}
plink --vcf babirusa_workshop_set.vcf.gz \
  --homozyg \
  --out babirusa_workshop_set_PLINK_A
```
Note that the backslash (`\`) allows you to continue writing the commands and not considering your command is done after you type it. This is just one way of doing it that makes it easy to read or follow up. You can also run it as `plink --vcf babirusa_workshop_set.vcf.gz --homozyg --out babirusa_workshop_set_PLINK_A` and it will run the same, just in one line.

The output files of the command will be a set of text files with the prefix declared in the `--out` argument ending with suffix "`.hom.*`". The "`.log`" file will contain all the options and the command.
  
```sh
ls -lh babirusa_workshop_set_PLINK_A.*
-rw-rw-r-- 1 sabhrina1 sabhrina1 687K Dec 13 08:03 babirusa_workshop_set_PLINK_A.hom
-rw-rw-r-- 1 sabhrina1 sabhrina1  827 Dec 13 08:03 babirusa_workshop_set_PLINK_A.hom.indiv
-rw-rw-r-- 1 sabhrina1 sabhrina1 161M Dec 13 08:03 babirusa_workshop_set_PLINK_A.hom.summary
-rw-rw-r-- 1 sabhrina1 sabhrina1 1.2K Dec 13 08:03 babirusa_workshop_set_PLINK_A.log
-rw-rw-r-- 1 sabhrina1 sabhrina1  170 Dec 13 08:03 babirusa_workshop_set_PLINK_A.nosex
```
  
Have a look on each of the files using `less` to see what each file contains. The `.hom` file contains detail of the ROH segments in the entire VCF. The `.hom.indiv` is the summary statistics of the ROH distribution per individual. Does all samples have ROH segment?

Read more about the output details in this [link](https://www.cog-genomics.org/plink/1.9/formats#hom).
  
The power of ROH is that the segment length is inversely correlated with the time of inbreeding event. The longer the segment is, the more recent is the inbreeding event. However, we cannot see this easily from the output file. We need to download the output files and plot the results. Before we do that, we need to prepare a R working directory in our local computer.
  
### Task 2: Plotting PLINK results in R
  
We will use R in R Studio to visualise our ROH results. To make our working directory, open R Studio, choose `File > New Project` and in the `New Project Wizard` choose `Create New R Directory > New Project`. Type in your directory name, such as "`Inbreeding_Analysis`" and choose your preferred path within your local computer. If you are happy with the directory name and location, click `Create Project`. A new R Studio session will appear.
  
To start our data visualization project, choose `File > New File > R Script`. Save the R Script from the start to avoid forgetting to save your commands. For example, save it as "`01_plot_ROH.R`" in our R project working directory.
  
Then, go to your command line interface and **go to the directory of your newly made R project**. Make a new directory called "`input`" and download your PLINK results there by doing `sftp` in `input`. Do not forget to also download the accompanying metadata.
```sh
mkdir input
cd input
sftp -i <path_to_identity_file> <username>@16.171.154.110
> get day_4_inbreeding/babirusa_workshop_set_PLINK_A.hom* .
> get day_4_inbreeding/babirusa_workshop_metadata.txt .
```
  
> Tips: You can check whether you are in the correct directory by looking at the content. If it contains a file with the directory name with .Rproj suffix, you are on the right place.
  
Afterwards, we go back to our R Studio and open our "`01_plot_ROH.R`" we just made a couple of paragraphs ago.

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
The argument `by=c("Sample"="IID")` allows you to merge similar columns. `Sample` is the column with the Sample IDs on `m` (the metadata) and `IID` is the column with Sample IDs on `p` (plink `.hom` result). Look at the result of the merge using `head()`. Is the merging successful? Pay attention in the observation count of the objects.

Note also that the number of observation has changed. This is because the samples in the metadata not all having ROH results while `.plink.hom` only gives ROH segment details results. For these samples, there are "`NA`" instead on each column. We want them to be "0" to reflect the lack of ROH later when plotting.
```
r<-replace(r, is.na(r), 0)
```

Now, we will look at how ROH segment distributed along chromosome 1. We can plot it using ggplot as follows:
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
  
Which region has the highest ROH? Could you say it from the plot?

> Tips: Along this exercise, feel free to change the arguments & exclude some lines to see what each part does. You can also do `?<function name>()` to see what a function does in R.
  
Our eyes can only make an estimate. Let's try to do a more certain visualization to do answer this question. First, we need to summarise the ROH length per sample. This can be done with `tidyverse` functionality as follows.
```{r group}
rSROH <- r %>%
  group_by(Sample) %>%
  summarise(SROH=sum(KB))
```
  
Open `rSROH` and see what is inside. You can use `cat()` or `head()`.
```
head(rSROH)
# A tibble: 6 × 2
  Sample    SROH
  <chr>    <dbl>
1 RD1    163709.
2 RD10   346205.
3 RD12        0 
4 RD16    42996.
5 RD17   243026.
6 RD2    179293.
```
This is the summary of the total amount of ROH segment we obtained for each sample. We have this summary on the `.hom.indiv.` Can you check whether we have the correct sum per sample? Use `head` or `less` command to peek on a handful of samples.
  
> Exercise 2
>
> Modify the above command by changing `sum` to `mean` to get the mean length of ROH segments detected. Then, check the `KBAVG` column in `.hom.indiv` to see if we have it correct.

To see which region has the most ROHs, we can merge `rSROH` with `r` and plot them as follows:
```{r}
r2<-left_join(r,rSROH,by=c("Sample"="Sample"))
  
ggplot(r2)+
  geom_boxplot(aes(x=Region,y=SROH/1e6,fill=Region))+
  labs(x="region",y="Total ROH amount (Mbp)")+
  theme_minimal()
```
  
Which population has the highest amount of genome segments indicative of recent inbreeding?
  
> Exercise 3
>
> The inbreeding coefficient (FROH) is the total sum of ROH segment divided by the entire genome length. As we map to the domestic pig genome and ran this on autosomal genome only, the total here is 2265774640 base pairs. Calculate the FROH and plot the result using a boxplot. The FROH should range from 0 to 1. Is it still looking the same with the SROH plot?
  
This plot, however, does not tell us how many are the short segments and long segments (remember the length of ROH is related to the time since inbreeding i.e. how old is the inbreeding). To do that, we will use `geom_freqpoly()` to show the distribution of ROH segment lengths in a sample.
```{r}
ggplot(r2)+
  geom_freqpoly(aes(x=KB,color=Region, group=Sample),
                  size=2, alpha=0.5)+
  xlim(1000,10000)+
  theme_minimal()
```

What could you conclude from this plot? What is the advantage of `geom_freqpoly()` compared to `geom_boxplot`?

The distribution of ROH segments has many features. Other than the total sum of ROH, its proportion across genome, and the mean length, there is also the number of ROH segments. We often called these summary statistics SROH, FROH, LROH, and NROH respectively. As ROH segments broken down with recombination across generations, segments from older inbreeding events will become shorter. As a result, a sample with a lot of ancient inbreeding event will have a lot of NROH, as we have shown in the lecture. Let's plot it.

First of all, we need to get the number of ROH segments per sample. We can use the same tidyverse functionality and do this:
```
rNROH <- p %>%
  group_by(IID) %>%
  summarise(NROH=n())
```

The `rNROH` object should be looking like this when seen with `head()`
```
# A tibble: 6 × 2
  IID    NROH
  <chr> <int>
1 RD1     124
2 RD10    251
3 RD16     34
4 RD17    170
5 RD2     135
6 RD20    152
```

Then, we merge `rNROH` with our previous data set
```
r3<-left_join(r2,rNROH,by=c("Sample"="IID"))
```

Then, we plot the spectrum of NROH and FROH using `geom_point()`
```
ggplot(r3)+
  geom_point(aes(x=FROH, y=NROH, color=Region), size=3)+
  labs(x="FROH", y="NROH")+
  theme_minimal()
```

What can you conclude from these plots?

> Exercise 4
> 
> For now we have ran PLINK with default options. These parameters values in the detault option were fine tuned for ROH analyses of human genomes. Now lets try to run PLINK with more sensible parameters for our species of interest (i.e. babirusa).
>
> The command is `plink --homozyg --homozyg-snp 20 --homozyg-kb 10 --homozyg-density 1000 --homozyg-window-snp 20 --homozyg-window-threshold 0.25  --vcf babirusa_workshop_set.vcf.gz --out babirusa_workshop_set_PLINK_B`.
>
> To understand what each argument does, read the documentation for `PLINK --homozyg` [here](https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo)
>
> After getting the output, download it via `sftp` and check the distribution of segment via plotting in R as we did with the previous `_PLINK_A` output. Are there difference in the results?

The parameters were coming from this paper, who used 60K SNPchip sequencing technology to get SNPs from domestic pig. Consequently, they were optimising PLINK to detect ROH on medium density SNPs (60000). In our VCF, we have more than 4 million SNPs representing the genome sequences.

<!---
I just added this, might not be a good idea? Just intense plotting session.
--->
  
### Task 3: Model-based method
  
As we just seen, changing parameters in PLINK can lead to widely different results. The most difficult issue for ROH detection is to distinguish true heterozygous sites from sequencing errors. Although sequencers today are doing a very good job they still make errors and these errors can creep in into ROH analyses in the form of false positive heterozygous sites which can result in surprious results (i.e. break ROH). In addition the amount of sequencing errors can vary depending on the sequencer, or even between runs of the same sequencer. 
  
Given we do expect some sequencing errors we need to allow for some "heterozygous" sites in ROH. Knowning how much we need to allow is an issue however. One approach to do this for example is to calculate "heterozygosity" on chromsome Y. Because chromsome Y is haploid we do not expect any heterozygous sites (except in the PAR region which is similar to a region of chromosome X). So the level of heterozygosity per Mb, computed on the Y chromosome, can be interepreted in as an expected number of false positive heterozygous sites per Mb. 
  
> Quick questions: do we expect the depth of coverage to be the same in autosomes and Y chromosome? How would that affect this analysis?
  
More complex, probabilistic model-based, methods have been designed to tackle this issue. This is the case or `ROHan` (https://github.com/grenaud/ROHan). ROHan works by calculating robust heterozygostiy estimate for specific regions of the the genome. It then uses a Hidden Markov Model (HMM) to identify ROHs, i.e. stretches of the genome that possess signficantly lower "heterozygosity" than expected as a result of: 1) heterozygosity that we can attribute to sequencing errors 2) or due to mutations that could have arisen since the individual inherited this segment of DNA.
  
`ROHan` probabilistic method is very long to run (~6 hours with 18 threads per sample!), so we have run it for you and we will focus on analysing the results instead. Here `ROHan` was run by setting the mutation rate in ROH parameter set to 0.0001. We will discuss more about this method during the tutorial's wrap up session.
  
The results of `ROHan` can be found in your `day_4_inbreeding` directory as you have made a symbolic link to all files in `/home/DATA/Day_4/Session_2`. These are the files starting with `RD` and ending with `.hmmrohl.gz`

> Exercise 5
>
> Have a look on one of the results using `zcat` and `head`.

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
  
Download these results into your R project directory where you plot your PLINK results. Then, read the files into R using `read.tables()`, add a column containing Sample ID, and concatenate all files into one file to make a similar dataframe as we have with PLINK where it contain all details of the ROH of all samples.
  
For example:
```{r}
s1<-read.table("input/RD1_mdup_1e-4.mid.hmmrohl.gz")
colnames(s1)<-c("rohID","chr","start","end","length","validatedSites")
s1$SampleID<-"RD1"
```
  
> Challenge: How can you make your task faster with R? Hint: make a list of the file paths.
  
If you do not feel challenged, feel free to repeat the commands manually or copy paste and run the following code (Warning: you need to install an optional package `janitor` and `gtools`).
```r
library(gtools)
library(janitor)
  
s<-list.files(path="input/",pattern="*.mid.hmmrohl.gz",full.names = T)  ### <---- the "path" needs to be adjusted to where you do the directory relative to the R script

head(s) ### <--- check first if your file path correct

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
  
The resulting file should look like this.
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
  
> Exercise 6
>
> Plot the segment distribution using `geom_segment()` as you have for the PLINK results. Which one is more similar to ROHan, PLINK_A or PLINK_B?
 
If you have done up to this point and there is still some time before we start the wrap up session, feel free to have another PLINK run with different parameters. For a start on sensible parameters and how to understands them, have a read on [this review paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6463-x) on how usually it is done with various domestic species.
