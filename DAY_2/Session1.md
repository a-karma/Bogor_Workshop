![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data pre-procesing - Session 1

Before we start, let's organize yesterday's material and prepare the directory structure for today.
Please connect to the server either via `putty` or via the `ssh` command then run:
```sh
mkdir day1; mkdir day2;
mv session2 ./day1/
mv project_bash ./day1/
cd day2
mkdir qc; mkdir bams; mkdir vcfs; mkdir fastqs; mkdir scripts; mkdir lists
touch what_i_did.txt
```
You should now be in your `~/day2/` directory, run `ls -lh` to ensure that all sub-directories have been created.

Let's now create a symbolik link to the directory containing the input data for the tutorial and activate the conda environment to access all software we will need:
```sh
ln -s /home/DATA/Day_2/ ~/day2/raw_data
conda activate Day_2
```
In your `raw_data` folder you should now see 8 files with the `.fastq` extension. These are the results of pair-end sequencing on Illumina HiSeq X platform of 4 babirusa individuals (one from each region of Sulawesi plus the Togean Islands as shown on the map). We are now going to familiarise with this bioinformatic file format and then evaluate the quality of these sequencing results.

![babirusa_map](../IM/babirusa_day2.png)

#### Fastq format
Genomic information from high-throughput sequencing are stored in text-based files called
FASTQ. These files comprise a series of entries containing not only the sequence but also
its quality score. Each entry consists of four lines as shown in the box below.

```sh
@HiSeq 4000:1:FCX:4:15:66:165 2:N:0:7
ATTTAGTACCATGACATGACACATACTACAATTGACGACATCAATCA
+
IGHFDEC@;;?=>B=?<;A:?@>9<>9756867544312*,*)'&)+
```

The first line (a.k.a header_line) starts with `@` and contains a series of characters which uniquely identify the read. 
When dealing with Illumina data, the identifier includes 7 fields encoding information about the sequencing process plus 4 fields about the read itself (see tab2.1).
Fields are separated by a `:` and the two groups are separated by a space.

Table 2.1: Fastq format – Fields in the Sequence ID line
| Field | example|
| ------ | ------ |
| Instrument used for sequencing|HiSeq 400|
| Run Number on that instrument|1|
| flow cell ID | FCX |
| lane number | 4 |
| tile number | 15 |
| X coordinate of cluster| 66 |
| Y coordinate of cluster| 165 |
| read number | 2 (2 nd read of the pair) |
| is filtered | N (Y did not pass, N otherwise) |
|control number| always zero on HiSeq X |
|sample number| 7 (the 7 th sample in the pool) |

The second line contains the actual sequence (in the example above, a short fragment of 47 nucleotides).

The third line always start with a `+` to improve the readability of the file by separating the
sequence line from its quality. The + symbol might be followed by the read identifier in
the first line (optional).

The forth line consists of a series of characters that encode the Phred-score for each
base in the sequence on the second line of the entry. Broadly speaking, Phred scores
are a measure of how confident we are in calling a given base.

Let P denote the probability of identifying the wrong nucleotide. Then, the Phred quality
scores Q are defined as:

Q = −10 log<sub>10</sub> P

Thus, if a base has an assigned Q-score of 20, it means that the chance that we have
called the wrong base are 1 in 100 i.e. the base call accuracy is 99%.

In the example above, the first base of the sequence is an A. The corresponding quality
of the base calling is encoded as `I` (first character on line 4). The symbol `I` corresponds
to a Q-score of 40 (see table 2.2) which means that the base call accuracy of the first
position in the sequence is 99.99%.

Table 2.2: Illumina Phred-score encoding
|Symbol|Q-score||Symbol|Q-score||Symbol|Q-score||Symbol|Q-score|
|----|----|----|----|----|----|----|----|----|----|----|
|!|0||,|11||7|22||B|33|
|”|1||-|12||8|23||C|34|
|#|2||.|13||9|24||D|35|
|$|3||/|14||:|25||E|36|
|%|4||0|15||;|26||F|37|
|&|5||1|16||<|27||G|38|
|’|6||2|17||=|28||H|39|
|(|7||3|18||>|29||I|40|
|)|8||4|19||?|30|
|*|9||5|20||@|31|
|+|10||6|21||A|32|


> Question 1: What's the accuracy of the last base call in the example above? 

#### Quality Control

Assessing the quality of sequencing results is a crucial step in genomic analysis. 
If your initial input are problematic, all downstream analysis might suffer from various biases and in general your inference will be less reliable.
In this tutorial we are going to use a ![FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) which is a program designed to spot potential problems in high througput sequencing datasets. 

The program will analyse the input file provided and it will produce a report (.html file) that can be visualised in a web-browser. Each tab in the FASTQC report will show different aspects of the quality of the sequencing but the `Per base sequence quality` is arguably the most important. 

> Question 2:
>
> In the figure below you can see two examples of very different quality control results. Can you guess which one is "the good run"?
 
![fastQC](../IM/fq_report.png)

Let's analyse our data and see what we get, shall we?

The command to run fastqc is quite straightforward:
```sh
fastqc -t 1 ./raw_data/sub_RD56_1.fastq -o qc
```
where `-t 1` corresponds to the number of threads used for the analysis, and `-o qc` corresponds to the output directory where the program will store the report.

In principle we could run this command for the other 7 files in our initial dataset by manually changing the file name every time. Given that this is quite tedious, we are going to make use of what we learnt yesterday and use a loop instead. First of all let's make a list of input files:

```sh
ls ~/day2/raw_data/*.fastq > ~/day2/lists/fastq_list.txt
```
You can inspect this list as usual using either the command `cat` (to print it on screen) or `less` (press q to exit)
Now that we have our list of files, we are going to run all the quality control analysis sequentially (one after the other) using our `while` loop:

```sh
while read -r line
do
fastqc -t 1 $line -o qc
done<./lists/fastq_list.txt
```
Once completed, use either the PSFTP app or the sftp command to transfer the .html files to your local computer and visualise them in your web-browser.
After establishing the sftp connection you can run:
```sh
cd ./day2/qc/
lcd path/to/download/folder
get *.html
```

> Question 3:
>
> How many reads we got for each individual?

#### Removing Adapters
After quality control, the next step in the pre-processing of NGS data consists in removing adapters. 
There are many software available that can perform this task, here we will focus on ![AdapterRemoval](https://adapterremoval.readthedocs.io/en/stable/). 
The main reason why we present this software in this workshop is because it not only search and remove adapters from high-throughput sequencing data 
but it can also perform the colllapsing of the two reads if necessary. 

Collapsing reads should be avoided when dealing with modern DNA given that the DNA molecules present in the library are longer than the read size (normally 150 bp for Illumina) hence the forward and reverse sequencing product do not overlap. Merging the overlapping region of the two reads is instead the standard procedure when analysing shorter DNA fragments which is always the case when dealing with degraded DNA such as DNA extracted from archaeological specimens or museum collections. 

In orther to remove adapters and low-quality bases at the termini of each read you can run the following command:

```sh
AdapterRemoval --file1 reads_1.fq --file2 reads_2.fq --basename output_paired --trimns --trimqualities
```
If working with ancient DNA (aDNA), you can add the `--collapse` option to merge the two reads and adjust the quality score accordingly.



