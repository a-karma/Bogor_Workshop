![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
please do not modify these first two lines of the .md file
use the same syntax to add pictures:
placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Detecting Recent Inbreeding using Runs of Homozygosity

### Introduction

In the last decade, ROH analyses have become the state-of-the-art method for inbreeding assessment. There are many ways to detect segments of ROH in a whole genome sequence. We can do it straight from BAM file (ROHan), but most do it from a VCF file. The principle is to slide a "window" across the VCF file and see whether there are subsequent homozygous genotypes stretches unbroken across the genome of each sample.

In this tutorial, we will use the most commonly used method: PLINK. PLINK is probably the most used program for analyzing SNP genotypes and ROH in human and animal populations.  In PLINK, the --homozyg function is used to perform ROH analyses and relies on several input settings.

### Task 0: Preparing your working directory

> Exercise 0.1
>
> Prepare your working directory for this session by following all the steps in Task 0 of Day 4 Session 1, but with data from /home/DATA/Day_4/Session_2

You should have a working directory in your home folder now named `day_4_inbreeding` containing a symbolic link of `babirusa_workshop_set.vcf.gz`.

### Task 1 Detecting ROH with PLINK

To detect homozygous segments in a VCF, we run a plink command as follows:
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

The detail of the ROH segments found can be found in the ".hom" file. (Read more about the output details in this [link](https://www.cog-genomics.org/plink/1.9/formats#hom)).

Let's have a look on the .hom.indiv to see the summary. Does all samples have ROH segment?

The power of ROH is that the segment length is inversely correlated with the time of inbreeding event. The longer the segment is, the more recent is the inbreeding event. However, we cannot see this easily from the output file. We need to download the output files and plot the results. Before we do that, we need to prepare a R working directory in our local computer.

### Task 2 Plotting PLINK results in R

We will use R in R Studio to visualise our ROH results. To make our working directory, open R Studio, choose File > New Project and in the `New Project Wizard` choose `Create New R Directory > New Project`. Type in your directory name, such as "Inbreeding_Analysis" and choose your preferred path within your local computer. If you are happy with the directory name and location, click `Create Project`. A new R Studio session will appear.

To start our data visualization project, choose `File > New File > R Script`. Save the R Script from the start to avoid forgetting to save your commands. For example, save it as "01_plot_PLINK.R" in our R project working directory.

Then, go to your command line interface and go to the directory of your newly made R project. Make a new directory called "input" and download your PLINK results there. Do not forget to also download the accompanying metadata you have sym-link-ed. For example:
```sh
mkdir input
cd input
sftp -i <path_to_identity_file> <username>@138.246.238.65
> get /home/<username>/day_4_inbreeding/babirusa_workshop_set_PLINK_A.hom* .
> get /home/<username>/day_4_inbreeding/babirusa_workshop_metadata.txt
```

> Tips: You can check whether you are on the correct directory by looking at the content. If it contains a file with the directory name with .Rproj suffix, you are on the right place.

Afterwards, we go back to our R Studio and open our "01_plot_PLINK.R".

To plot the results, we need to first read the PLINK results and the accompanying metadata.
```{r readHom}
library(dplyr)
p<-read.table("input/babirusa_workshop_set_PLINK_A.hom", header=T)
m<-read.table("input/babirusa_workshop_metadata.txt",header=T)
```

To check whether our dataset were properly read into the assigned object `p` and `m`, we use the R command `head()`.
```{r}m<-read.table("babirusa_workshop_metadata.txt",header=T)
head(p)
head(m)
```

We would like to plot while including the metadata. To do this, we will use `left_join()` function from the package `dplyr` to join the two.
```
r<-left_join(m,p,by=c("Sample"="IID"))
```
The argument `by=c("Sample"="IID")` allows you to merge similar columns. Look at the result using `head()` is the merging successful? Pay attention in the observation count of the objects.

Now, we will look at how ROH segment distributed along chromosome 1? We can plot it using ggplot as follows:
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

Which region has the highest ROH? To know that, we need to summarise the ROH length per sample. This can be done with tidyverse functionality as follows.
```{r group}
rSROH<-r %>% group_by(Sample) %>% summarise(SROH=sum(KB))
```
Open rSROH and see what is inside. This is the summary of the total ROH segment we obtained.

> Exercise 2.1
>
> We have this summary on the .hom.indiv. Can you check whether we have the correct sum per sample?
> Repeat the above command to get the mean length and also check.

To see which region has the most ROHs, we can merge `rSROH` with `r` and plot them as follows:
```{r}
r2<-left_join(r,rSROH,by=c("Sample"="Sample"))

ggplot(r2)+
  geom_boxplot(aes(x=Region,y=SROH/1e6,color=Region))+
  labs(x="region",y="Total ROH amount (Mbp)")+
  theme_minimal()
```

Which population has the highest amount of inbreeding?

> Exercise 2.2.
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

Let's go back to the server and do Exercise 2.4.

> Exercise 2.4.
>
> We have run the default PLINK option that was based on human genome. Try to run PLINK with minimum SNP count 20, minimal ROH segment length 10 kb, 1 SNP per 1 Mbp, and scanning window consists of 20 SNP with 0.25 portion of the overlapping windows must be called homozygous to define any given SNP as 'in a homozygous segment'.
> Hint: Read the documentation for PLINK --homozyg here https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo
>
> Compare the result of Exercise 2.4 with the result from the default. Which has more ROHs per sample? Which set of parameters do you trust more?
