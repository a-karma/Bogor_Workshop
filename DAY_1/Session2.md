![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
## Day 1 - Basic concepts of command line programming - Session 2

### 1. Shell Scripting
In the previous session we have seen some


### 2. Working with bioinformatic softwares using conda
In Session 1 we have seen three examples of text file that are commonly used in bioinformatics:

- The `fasta` format to store DNA sequence information
- The `GTF` (General Transfer Format) developed for the Ensembl genome browser
- The `BED` (Browser Extensible Data) format developed at UCSC for the Genome Browser tool

The `GTF` and the `BED` format are both TAB-separated files used to store genomic regions as coordinates along with their associated annotations. 
Although in principles we could manually edit these files using standard text editors this becomes very unpractical when dealing with very large files.
A more efficient option would be combinig command line tools (like `sed`, `awk`, or `grep`) but performing complex tasks using only these tools is not straightforward and often require a deep understanding of programming. Luckly for us, bioinformaticians have created various software specifically designed to manipulate these file formats. 

Let's have a look at `bedtools` a powerful toolset for genome arithmetic. 

In your terminal, please type: 

```sh
bedtools --help
```

Looks like the program is not installed :( 

To protect the integrity of the file system on a server, normal user do not have permissions to directly install softwares. Moreover, almost any bioinformatic tool will rely on specific libraries or package versions that might create conflicts or even impeed the functionality of other programs. To circumvent these issues, we need to make sure that the software we need are installed in a "confined space" (a.k.a an environment) containing all the necessary dependencies that can become accessible only when we need it. This can be esily implemented using `conda` which is a package, dependency, and environment manager for any programming language. You can read more about conda [here](https://docs.conda.io/en/latest/).

We have already installed conda on our server and we have created a different environment for each day of the workshop. In order to have access to conda you need to run:

```sh
source /home/anaconda3/bin/activate
conda init
```
You need to run this command only once and you should see a change in your prompt: the word `(base)` appears on the left.
This is telling us that you are now in the conda base environment which represent the default space. 
To activate the environment for this session, please run:

```sh
conda activate Day_1
```

> Question: Have a look at your prompt again, what do you see? 


Now that bedtools is accessible let's see what it can do.

Suppose you are interested in analysing only neutral evolving sites. Thus, you may want to remove sites that are likely to be under selective pressures. 
As a first approximation you can start to analyse SNPs that do not fall inside CDS. You can just run the following commands:

```sh 
bedtools intersect -a snps_panel.bed -b genes_chr30.gtf -v > snp_filtered.bed
```
The `-v` flag tells `intersect` to report all lines in the target file (specified using the `-a` flag) that do not overlap with the genomic intervals listed in the bed file (-b flag).

There is a lot more that you can do with genome arithmetic. Let’s picture a more complex
scenario: suppose you are interested in studying the promoter regions of various genes
on this chromosome. You have a fasta file with the entire chromosome sequence (see
the ptw_ch30.fa file) and you would like to examine 10 kb upstream the starting codon
of each gene. How can you get that information?

First of all you need to get the coordinates of the starting codons. This is very easy using grep:
```sh
grep 'start_codon' genes_chr30.gtf > CDS_start.gtf
```

Now you need to modify this file using the function `flank` implemented in bedtools which will flanking intervals for each region in a BED/GFF/VCF file.
This function requires also a genome file defining the length of each chromosome, so let’s create this file first.
```sh
echo -e "chr30\t150000000" > ch30_length.bed
```

The above command `echo` would normally output on screen whatever string you type after it. 
The `-e` option tells the software to enable the interpretation of backslash escapes and `\t` stands for TAB.
Now we can run:

```sh
bedtools flank -i CDS_start.gtf -g ch30_length.bed -l 10000 -r 0 > ch30_promoters.gtf
```
The command above creates a new 
Finally we can use the `getfasta` function in bedtools to extract the sequences of the promoter regions:

```sh
bedtools getfasta -fi ptw_ch30.fa -bed promoter.bed -fo ptw_prom_sequences
```

In the command above the `-fi` flag stands for file input, the `-bed` indicates the coordinate file while the `-fo` option stands for file output. 
You can examine the first output line using: 
```sh
cat ptw_prom_sequences | head -1
```
> Exercise 2:
>
> Try to combine the intersect and flank functions in order to filter the `snp_ch30.bed` file
> by excluding CDS and all reagions that are 5 kb form the starting and the stop codon of each CDS.

### 3. Transferring files
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
