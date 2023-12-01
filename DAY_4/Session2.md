![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
please do not modify these first two lines of the .md file
use the same syntax to add pictures:
placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Introduction

There are many ways to detect segments of ROH in a whole genome sequence. We can do it straight from BAM file (ROHan), but most do it from a VCF file. The principle is to slide a "window" across the VCF file and see whether there are subsequent homozygous genotypes stretches unbroken across the genome of each sample.

In this tutorial, we will use the most commonly used method: PLINK in VCF and BCFtools ROH. While the first use observational method relying on counting homozygous and heterozygous genotypes across certain base pair range, the second detects regions of autozygosity in sequencing data using a hidden Markov model. 

### Task 1 Detecting ROH with PLINK

To detect homozygous segments in a VCF, we run a plink command as follows.
```{bash plink, eval=FALSE}
plink --vcf babi_set4_qualFilt_geno0_maf0.05_LD0.1.vcf --homozyg --out babi_set4_qualFilt_geno0_maf0.05_LD0.1_Bosse2012 --homozyg-density 1000 --homozyg-kb 10 --homozyg-snp 20 --homozyg-window-snp 20 --homozyg-window-threshold 0.25
```

The output files of the command will be a set of text files with the prefix declared in the --out argument ending with suffix ".hom.*". The detail of the ROH segments can be found in the ".hom" file.

```{r readHom}
library(dplyr)
pl<-read.table("input/babi_set4_qualFilt_geno0_maf0.05_LD0.1_Bosse2012.hom", header=T)
head(pl)

m<-read.table("babirusa_workshop_metadata.txt",header=T)
head(m)

r<-inner_join(m,pl,by=c("Sample"="IID"))
glimpse(r)
```

How does the ROH segment distributed along chromosome 1? We can plot it using ggplot as follows:
```{r plotSeg}
library(ggplot2)
r %>%
  filter(CHR==1) %>%
  ggplot()+
  geom_segment(aes(x=POS1, xend=POS2,
                   y=Sample, yend=Sample,
                   group=Sample, color=Region), size=3)+
  labs(aes(x="Chromosome 1 (Mbp)",y="Sample ID"))+
  theme_minimal()
```

We can plot the summary statistics of ROH segments as follows.
```{r readHomInd}
pl2<-read.table("input/babi_set4_qualFilt_geno0_maf0.05_LD0.1_Bosse2012.hom.indiv", header=T)
r2<-inner_join(m,pl2,by=c("Sample"="IID"))
glimpse(r2)
```

```{r}
ggplot(r2)+
  geom_boxplot(aes(x=Region, y=KB,color=Region))+
  labs(aes(x="Region",y="Total ROH segment"))+
  theme_minimal()
```

Question: Which population has the highest amount of ROH segment?

### Task 1b (Optional)

Now try the default option. Do not forget to make the --out prefix different from the previous command.
```{bash plink_def, eval=FALSE}
plink --vcf input --homozyg --out output_def
```

Question: Have a look on the files. What is different and why is it so? You can have a look on the full list of default setting in this link: https://www.cog-genomics.org/plink/1.9/ibd#homozyg

### Task 2 Detecting ROH with Bcftools ROH

To use `bcftools roh` we run the following command:
```{bash bcftools, eval=FALSE}
bcftools roh --AF-dflt 0.4 -G 30 babi_set4_qualFilt_geno0_maf0.05_LD0.1.vcf -o babi_set4_qualFilt_geno0_maf0.05_LD0.1_ROH
```

You will wait for less than 5 minutes and obtained the following message from the software:
```{bash bcftools_out, eval=FALSE}
Number of target samples: 18
Number of --estimate-AF samples: 0
Number of sites in the buffer/overlap: unlimited
Number of lines total/processed: 200976/200976
Number of lines filtered/no AF/no alt/multiallelic/dup: 0/0/0/0/0
```

The resulting output contains two parts. The first one, RG, stands for "ROH over a region". The "ST" stands for "state" which gives more details on whether a SNP is on Hardy Weinberg ("HW") or Autozygous ("AZ").

We need to first separate the output as follows:
```{bash bcftools_out_RG, eval=FALSE}
grep "RG" babi_set4_qualFilt_geno0_maf0.05_LD0.1_ROH > babi_set4_qualFilt_geno0_maf0.05_LD0.1_ROH.RG 
```

Then, we downloaded the results to our local computer for plotting.
```{bash file_transfer, eval=F}
sftp -i <path_to_your_identity_file> <your_username>@138.246.238.65
```

To plot the ROH results, we first need to read the file into R
```{r readFile}
rg<-read.table("input/babi_set4_qualFilt_geno0_maf0.05_LD0.1_ROH.RG")
colnames(rg)<-c("RG", "SampleID", "chr", "start", "end", "length", "markers", "quality")

head(rg)
```

Note that the Sample ID is a double. We can first remove the extra sample ID before adding the metadata.
```{r}
rg$SampleID<-gsub("_RD[0-9]*","",rg$SampleID)
head(rg)

br<-left_join(m,rg,by=c("Sample"="SampleID"))
head(br)
```

There is an info of marker density per detected ROH segment length here. Let's use that to filter for ROH with low SNP support.
```{r density}
br<-br%>%mutate(density=length/markers) 

ggplot(br)+
  geom_histogram(aes(x=density))+
  theme_minimal()

ggplot(br)+
  geom_histogram(aes(x=quality))+
  theme_minimal()
```

Based on the histogram, filtering on quality of 20 and density of 10000 might be reasonable.

Let's use it to detect distribution of ROH segment.
```{r rgPlot_segment}
br %>%
  filter(density>10000) %>%
  filter(quality>20) %>%
  ggplot()+
  geom_segment(aes(x=start, xend=end,
                   y=Sample, yend=Sample, color=Region), size=3)+
  theme_minimal()
```

```{r}
br %>%
  filter(quality>20) %>%
  filter(density>10000) %>%
  ggplot()+
  geom_boxplot(aes(x=Region, y=length,color=Region))+
  labs(aes(x="Region",y="Total ROH segment"))+
  theme_minimal()
```

Question: Why do you think bcftools and PLINK differs and which one do you trust more?
