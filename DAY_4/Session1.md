![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

<!---
[comment]: # please do not modify these first two lines of the .md file
[comment]: # use the same syntax to add pictures:
[comment]: # placeholder name within square brackets and ../IM/file_name.png within parentheses
--->

## Inferring Ancient Demography with PSMC

### Introduction

The Pairwise Sequentially Markovian Coalescent (PSMC) model uses information in the complete diploid sequence of a single individual to infer the history of population size changes. Originally published in 2011 (Li and Durbin 2011), it has become a very popular tool in the world of genomics. In this tutorial, we walk through the steps to generate the necessary input data for PSMC and run it on our babirusa data set to see how well we can estimate effective population size. 

We will use the `psmc` tool to estimate the effective population size for our data. Instead of directly using bam files, we will use a variant calling format (vcf) file as the starting point. 

### Task 0: Preparing your working directory

Before starting this first session, activate the conda environment containing all the software we are going to use today as follows:

```sh
conda activate Day_4
```

In this exercise, you will need the following input files which you can find in `/home/DATA/Day_4`:
```sh
babirusa_workshop_set.vcf.gz
babirusa_workshop_set.vcf.gz.tbi # this is just an index file to make bcftools works easily
chr_autosomes.txt
```

Make your own directory for this project and use a symbolic link for easy access to your working directory:
```sh
mkdir day4_psmc_tutorial
cd day4_psmc_tutorial
ln -s /home/DATA/Day_4/babirusa_workshop_set.vcf.gz* .
ln -s /home/DATA/Day_4/chr_autosomes.txt .
```

Now we are good to go.

### Task 1: Preparing input file for PSMC

PSMC works on a psmcfa file of a single diploid sample. A psmcfa file is a fasta representation of the genome (or part of a genome), where each entry (letter) represents if there exists a heterozygous variant in a fixed size window. Let us try and understand this by converting our vcf file to a psmcfa file. 

As we have a multisample VCF, we need to first extract a single diploid sample from this VCF. In this example, we will work with sample RD44 and RD71. The command for sample RD44 has been given; the command for sample RD71 is left for participants for their exercise.
```sh
bcftools view -s RD44 babirusa_workshop_set.vcf.gz > RD44.vcf
bgzip RD44.vcf
tabix RD44.vcf.gz
```
Do not forget to compress the resulting single sample VCF and index it with `tabix`

> Food for thought: Can you guess why RD44 is chosen? Hint: Have a look on the metadata.

Then, we create a consensus sequence using `samtools faidx` and `bcftools consensus`. We will use the autosomes list on `chr_autosomes.txt` to keep our analysis within the autosomal chromosomes.
```sh
REF=/home/DATA/Day_3_b/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
samtools faidx $REF -r chr_autosomes.txt | bcftools consensus RD44.vcf.gz > RD44.fa
```

> Food for thought: Why don't we use mitochondrial and sex chromosome?

After getting our consensus fasta, we construct our psmcfa file using `fq2psmcfa`.
```sh
fq2psmcfa RD44.fa > RD44.psmcfa
```

<details open>
<summary>Have a look on the result using `head` command. What do you see?</summary>
<br>
It appears to be a regular fasta file, with only "T" and "K". Here a "T" represents a 100 bp window without any heterozygous sites in it, whereas a "K" represents a 100 bp window with at least 1 heterozygous site in it.
</details>

By this stage, you should already have RD44.psmcfa and RD71.psmcfa to run.

### Task 2: Running PSMC

Now it is time for us to run our first psmc analysis. Let us first look at the options for running PSMC by running `psmc` on the command line.

```sh
psmc
```

In particular, note the -p option, which allows us to specify a pattern of parameters. This parameters allow us to divide up time into discrete bins, and the pattern we specify in -p option allows us to estimate the same parameter for multiple consecutive time bins. 

For example - the default pattern "4+5*3+4" splits time up into 23 bins (4+15+4), where only 1 parameter is estimated for the first 4 bins, then 1 parameter is estimmated for the next 5 triples of bins, and 1 parameter for the last 4 bins. So the total number of paramters we are estimating is 7 (1+5+1). Note this also in the plots for the effective population size, when we plot them.

Let us now run psmc for the first time - we will just use the default values for the options, while still explicitly specifying the parameter pattern. This pattern is quite coarse, but we will try and estimate it again with finer time bins in the next section.

```sh
psmc -p "4+5*3+4" -o RD44_coarsePattern.psmc RD44.psmcfa &
psmc -p "4+5*3+4" -o RD71_coarsePattern.psmc RD71.psmcfa &
```

As we have done this for whole autosomal genomes, this will take ~20 minutes. While we wait for it to finish, let's try have a look on only one chromosome to have faster result.

### Task 3: Running PSMC in chromosome 1

To run PSMC on one chromosome only, we need to prepare a psmcfa for only one chromosome.

<details open>
<summary>Can you guess how it is done?</summary>
<br>
  See the steps on Task 1. We need to first single out chromosome one as we make the consensus fasta and made a .psmcfa.
  ```sh
  REF=/home/DATA/Day_3_b/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
  samtools faidx $REF 1 | bcftools consensus RD44.vcf.gz > RD44_chr1.fa
  fq2psmcfa RD44_chr1.fa > RD44_chr1.psmcfa
  ```
</details>

After having the psmcfa of chr1 from both samples, run the same command as Task 1:
```sh
psmc -p "4+5*3+4" -o RD44_1_coarsePattern.psmc RD44_chr1.psmcfa &
psmc -p "4+5*3+4" -o RD71_1_coarsePattern.psmc RD71_chr1.psmcfa &
```

This should run in only within a minute. Let's have a look on one of the output files:
```sh
tail -34 RD44_1_coarsePattern.psmc
```

> Exercise 3.1:
>
> Modify the -p option to have finer resolution. Use "4 + 10*2 + 4 + 6" and “4 + 25*2 + 4 + 6”.
> Name the resulting files differently. For example, adding suffix "fine_a.psmc" and "fine_b.psmc".

### Task 4: Plotting PSMC

Let us now cat the psmc outputs for all 2 samples, and plot them into a pdf.
```sh
cat RD44_1_coarsePattern.psmc \
    RD71_1_coarsePattern.psmc > combined_coarsePattern.psmc 
```

To simplify the psmc results, we need to run the .psmc file with `psmc_plot.pl`. First let have a look on this script to see what it does.
```sh
psmc_plot.pl
```

An important thing here is to make sure we use the correct mutation rate and generation time. Because there is no mutation rate estimates yet for babirusa, we used domestic pig mutation rate here and average generation time we found on babirusa in captivity.
```sh
psmc_plot.pl -u 1.5e-09 -g 3 -s 100 -Y 1 -m 5 -n 30 -p -M "SE, TO" babirusa_chr1_coarse combined_coarsePattern.psmc 
```
Note that "SE" and "TO" stands for the different population assignment of RD44 and RD71 respectively.

> Exercise 4.1:
>
> Re-run the entire command for Task 4 with the results of Exercise 3.1.
> Do not forget to choose different prefix name for files from different parameters!
> For example, babirusa_chr1_fine_a and babirusa_chr1_fine_b

By this stage, you should have three pdf files, each coming from different -p option. To see how they look like, we need to download the output into our local machine using sftp:
```sh
sftp -i <path_to_identity_file> username@138.246.238.65
> get babirusa_* .
```

> Exercise 4.2:
>
> Re-run commands from Task 3 and 4 to the resulting whole genome psmc you did in Task 2. Compare the result of this coarse pattern. How are they different?

### Challenge (Optional)

If you have done the above exercises, feel free to try working on other samples by modifying all the commands from Task 2 to Task 4 with different sample names. Choose samples of the highest coverage and see how they compare. Work only on one chromosome first to get faster results.
