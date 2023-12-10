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
  
  > Exercise 1
  >
  > Prepare your working directory for this session by following all the steps in Task 0 of Day 4 Session 1.
  >
  > This time but with data from /dev/workshop_DATA/Day_4/Session_2
  
  You should have a working directory in your `/de/workshop_participants/your_user_ID/` folder now named `day_4` with a subdirectory called `inbreeding` containing a symbolic link of `babirusa_workshop_set.vcf.gz`, `babirusa_workshop_metadata.txt`, and 18 files of a ROH run results ending with `.hmmrohl.gz`. We will work with these results a bit later in the session.
  
  ### Task 1: Detecting ROH with PLINK
  
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
  
  ### Task 2: Plotting PLINK results in R
  
  We will use R in R Studio to visualise our ROH results. To make our working directory, open R Studio, choose File > New Project and in the `New Project Wizard` choose `Create New R Directory > New Project`. Type in your directory name, such as "Inbreeding_Analysis" and choose your preferred path within your local computer. If you are happy with the directory name and location, click `Create Project`. A new R Studio session will appear.
  
  To start our data visualization project, choose `File > New File > R Script`. Save the R Script from the start to avoid forgetting to save your commands. For example, save it as "01_plot_PLINK.R" in our R project working directory.
  
  Then, go to your command line interface and go to the directory of your newly made R project. Make a new directory called "input" and download your PLINK results there. Do not forget to also download the accompanying metadata you have sym-link-ed. For example:
  ```sh
  mkdir input
  cd input
  sftp -i <path_to_identity_file> <username>@138.246.238.65
  > get babirusa_workshop_set_PLINK_A.hom* .
  > get babirusa_workshop_metadata.txt
  ```
  
  > Quick question: how do you download the entire directory in one command?
  
  > Tips: You can check whether you are in the correct directory by looking at the content. If it contains a file with the directory name with .Rproj suffix, you are on the right place.
  
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
  
  > Exercise 2
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
  
  > Exercise 3
  >
  > The inbreeding coefficient (FROH) is the total sum of ROH segment divided by the entire genome length. As we map to the domestic pig genome and ran this on autosomal genome only, the total here is 2,265,774,640 base pairs. Calculate the FROH and plot the result. Is it still looking the same with the SROH plot?
  
  This plot, however, does not tell us how many are the short segments and long segments (remember the length of ROH is related to the time since inbreeding i.e. how old is the inbreeding). To do that, we will use `geom_freqpoly()` to show the distribution of ROH segment lengths in a sample.
  ```{r}
  ggplot(r2)+
    geom_freqpoly(aes(x=KB,color=Region, group=Sample),
                  size=2, alpha=0.5)+
    xlim(1000,10000)+
    theme_minimal()
  ```

  <!---
Let's go back to the server and do Exercise 2.4.
  
  > Exercise 2.3.
  >
  > For now we have ran PLINK with default options. These parameters values in the detault option were fine tuned for ROH analyses of human genomes. Now lets try to run PLINK with more sensible parameters for our species of interest (i.e. babirusa). Lets set minimum SNP count to 20, minimal ROH segment length to 10 kb, a maximum of 1 SNP per Mbp. ####LF COMMENT: NOT SURE WHAT THIS SENTENCE MEANS HERE#### We set the size of windows to 20 SNPs with 0.25 portion of the overlapping windows must be called homozygous to define any given SNP as 'in a homozygous segment'.  Name the output files of this new run with `_PLINK_B` suffix.
  > 
  > Hint: Read the documentation for PLINK --homozyg here https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo
  >
  > Compare the result of Exercise 2.4 with the result from the default. Which analysis detect the most ROHs per sample? Which set of parameters do you trust more?
  --->
  
  <!---
  I just added this, might not be a good idea? Just intense plotting session.
  --->
  
  ### Task 3: Model-based method
  
  As we just seen, changing parameters in PLINK can lead to widely different results. The most difficult issue for ROH detection is to distinguish true heterozygous sites from sequencing errors. Although sequencers today are doing a very good job they still make errors and these errors can creep in into ROH analyses in the form of false positive heterozygous sites which can result in surprious results (i.e. break ROH). In addition the amount of sequencing errors can vary depending on the sequencer, or even between runs of the same sequencer. 
  
  Given we do expect some sequencing errors we need to allow for some "heterozygous" sites in ROH. Knowning how much we need to allow is an issue however. One approach to do this for example is to calculate "heterozygosity" on chromsome Y. Because chromsome Y is haploid we do not expect any heterozygous sites (except in the PAR region which is similar to a region of chromosome X). So the level of heterozygosity per Mb, computed on the Y chromosome, can be interepreted in as an expected number of false positive heterozygous sites per Mb. 
  
  > Quick questions: do we expect the depth of coverage to be the same in autosomes and Y chromosome? How would that affect this analysis?
  
  More complex, probabilistic model-based, methods have been designed to tackle this issue. This is the case or `ROHan` (https://github.com/grenaud/ROHan). ROHan works by calculating robust heterozygostiy estimate for specific regions of the the genome. It then uses a Hidden Markov Model (HMM) to identify ROHs, i.e. stretches of the genome that possess signficantly lower "heterozygosity" than expected as a result of: 1) heterozygosity that we can attribute to sequencing errors 2) or due to mutations that could have arisen since the individual inherited this segment of DNA.
  
  `ROHan` probabilistic method is very long to run (~6 hours with 18 threads per sample!), so we have run it for you and we will focus on analysing the results instead. Here `ROHan` was run by setting the mutation rate in ROH parameter set to 0.0001. We will discuss more about this method during the tutorial's wrap up session.
  
  The results of `ROHan` can be found in your `day_4_inbreeding` directory. Let's have a look on one of the results using `zcat` and `head`.
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
  
  If you do not feel challenged, feel free to repeat the commands manually or run the following code (Warning: you need to install an optional package `janitor` and `gtools`).
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
  
  The resulting file should look like this:
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
  
  > Exercise 4
  >
  > Plot the segment distribution using `geom_segment()`, a boxplot, and the frequency of segment length classes as you have for the PLINK results. How similar it is with the PLINK results?
 
If you have done up to this point and there is still some time before the wrap up session, feel free to have another PLINK run with different parameters. Read the documentation for PLINK `--homozyg` [here](https://zzz.bwh.harvard.edu/plink/ibdibs.shtml#homo) to have an idea of how the parameters can be tuned in. For a start, have a look on [this review paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6463-x) on how usually it is done with various domestic species.
