![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data procesing - Session 2
In this session you will learn how to align reads to a reference genome and how to manipulate bam files.

#### The SAM and BAM format
SAM, or Sequence Alignment Map, is a standard file format for storing and exchanging alignment data generated from high-throughput sequencing technologies. It is a tab-delimited text format that consists of a header section and an alignment section. The header section provides metadata about the alignment file, while the alignment section contains detailed information about each individual read's alignment to a reference genome.

The BAM format, or Binary Alignment Map, is a compressed binary version of a SAM file. It is more compact and efficient for storing large datasets, making it the preferred format for most downstream analysis.

Each SAM alignment line contains eleven mandatory fields, each representing a specific aspect of the read's alignment. These fields include information such as the read name, chromosome, strand, mapping position, alignment score, and mismatches. Additional optional fields can also be included to store more detailed information or aligner-specific annotations.

The SAM format is a widely used standard in genomics research, allowing for efficient communication and analysis of alignment data across different software tools and platforms.

In your raw_data directory you should see a folder called CHR_10. 
This contains a five bam files that we have generated for you using the same procedure that you are learning today. 
Let’s have a look at a bam file and look at it with the widely used tool called `samtools`.

First try to type:
```sh
less ~/day2/raw_data/CHR_10/RD10_chr10.bam
```
and as usual type Q to close it

This is how a binary file looks like, not very user friendly or human readable right?

Let’s find a better way to visualize it. We should start by examining the header section:

```sh
samtools view -H ~/day2/raw_data/CHR_10/RD10_chr10.bam | head -20
```
The `view` command of samtools allows you to read, print, and convert SAM/BAM/CRAM files. The `-H` flags outputs only the header.
As you can see, there are 19 lines starting with `@SQ`. Each of these lines corresponds to a chromosome of the reference genome and the value you read after `LN:` corresponds to the chromosome length in base pairs (bp).  

Now let’s have a look at an alignment line. We can use again the view command:
```sh
samtools view ~/day2/raw_data/CHR_10/RD10_chr10.bam | head -1
```
As you can see there are two fields occupying the largest portion of the line: one is the sequence (10 th field) and the other one is the Phred-scaled base quality (11 th field).
Those correspond respectively to the second and 4th line of the fastq file used to generate this bam. This is definitely better than using less but what if we wanted to visualize a specific region of this chromosome?

You may have noticed that in the CHR_10 folder there are also 5 files with the .bai extension. 
Those are the indexes of the bams and they allow fast access to specific regions without going through the whole alignment.
thanks to the indexing procedure we can use another samtools command and look at any given portion of chromosome 10:
```sh
samtools tview ~/day2/raw_data/CHR_10/RD10_chr10.bam -p 10:21100
```

We are looking at the alignment on chromosome 10 starting at the position 21100 from
the beginning of the chromosome (the `-p` flag stands for position). The second line is full
of Ns because we haven’t provide any reference genome in a fasta format to this samtools commad. 
The third line is the consensus sequence while each of the following line represent a specific read that was mapped to that region of the reference genome. 
As you can see, reads are displayed in different colours which represent quality scores. Just type `?` to look at the visualization options. Press for example `n` or `c` and see what happens.

If you are eager to know more about BAM files and samtools please try this ![tutorial](https://sandbox.bio/tutorials?id=samtools-intro) after the course. 

Now that we have familiarised ourselves with this file format, let's try to generate some new bams!

#### Aligning reads
Aligning short-read sequences is the foundational step to most genomic and transcriptomic analyses. 
There are several tools to perform this task such as Bowtie2, BWA, HISAT2, MUMmer4, STAR, and TopHat2 just to name a few.
Here we are going to focus on the Burrows Wheeler Aligner (BWA) but the procedure can be generalised to any aligner with minor modifications.

The first step for using BWA consist in makeing an index of the reference genome in fasta format. This can be done via the `bwa index` command:

```sh
bwa index [-a bwtsw|is] input_reference.fasta index_prefix
```

If you look into your raw_data//SUS_REF/ folder you will see that contains the sus scrofa reference genome (.fa) plus a few other files with various extensions. Those are exactly the index files produce by the above command.

The next step consists in choosing the aligner algorithm and map the reads to the reference genome. Here we will use `bwa mem` which performs local alignment and produces alignments for different part of the query sequence. The basic usage of bwa mem is:

```sh
bwa mem reference_index_prefix reads_pair_1.fastq reads_pair_2.fastq <bwa options>
```
The full list of `bwa mem` option can be found in the manual ![here](https://bio-bwa.sourceforge.net/bwa.shtml).

In our case we will use the `-t` option to specify the number of threads used for the computation and the `-R` flag to specify the read group information.
The output of bwa mem is a SAM file, in order to reduce disk usage we would like to produce directly its compressed binary version (i.e. BAM). We can easily achieve this by piping the output of bwa into samtools view. Here is an example of the full command that we are going to run for all four babirusa individuals that we have sequenced:

```sh
bwa mem raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa fastqs/read1_file fastqs/read2_file -t 1 -R read_group_info | samtools view -Shu - > bams/bam_file_name.bam

```
In the command above, the `read_group_info` correspond to a string starting with `@RG`. The string (enclosed in quotes) should contain the following tab-separated fields:

- ID = read group unique IDentifier. 
- PU = Platform Unit i.e. the flow cell and lane barcode (optional)
- SM = sample name
- PL = PLatform technology
- LB = LiBrary identifier.

Recording this information in our bams is crucial because it will allow us to identify and mitigate batch effects. Imagine a scenario in which each sample has been sequenced twice on two different flow cells. During the library preparation for the second run, two labels have been swapped but we have discovered the errror only after the sequencing data has been processed. do we need to throw away all our work and start again from scratch? Absolutely not! We can easily extract the wrongly assigned reads using the read group info and save a lot of computing time. 

Let's now see how to extract the required read group info from a fastq header:
```sh
head -1 ~/day2/raw_data/sub_RD56_1.fastq
```
The output on your screen should be:
```sh
@A00155:379:HMTKMDSXY:3:2406:8594:22608 1:N:0:GTTCCAAT+GCAGAATT
```
The relevant info for the ID/PU field are the flow cell id (HMTKMDSXY in our case) and the lane number (3). For the LB field it is common practice to use the the last portion of the fastq header containing the external indexes of the library (GTTCCAAT+GCAGAATT). Therefore for this sample our read_group_info should look like this:

```sh
"@RG\t\tID:HMTKMDSXY.3\tPL:illumina\tSM:RD56\tLB:GTTCCAAT+GCAGAATT"
```

> `Exercise 1`
>
> prepare a file with the read group information for all our samples (one sample per line). store this file in the `lists` directory as `rg_info.txt`

In order to perform the alignment step for all our samples we are going to use the same loop structure you have seen in session 1. Thus, we need to prepare the full argument list using `paste`.

```sh
ls ~/day2/fastqs/*pair1.truncated > ./lists/r1.txt
ls ~/day2/fastqs/*pair2.truncated > ./lists/r2.txt
paste ./lists/r1.txt ./lists/r2.txt ./lists/rg_info.txt > ./lists/bwa_full_arg_list.txt
```
After this we need to create our script for bwa and trannsform into an executable:

```sh
touch ./scripts/bwa_aligner.sh
chmod 770 ./scripts/bwa_aligner.sh
```
and edit it with nano. Here's the content of the script:

```sh
#!/usr/bin/bash
INPUT1=$1
INPUT2=$2
INPUT3=$3
OUTPUT=$(echo `basename ${INPUT1}` | sed 's/.pair1.truncated//')
bwa mem raw_data/SUS_REF/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa ${INPUT1} ${INPUT2} -t 1 -R ${INPUT3} | samtools view -Shu - > bams/${OUTPUT}
```

Let's now run our while loop (this will take a few minutes):

```sh
while read -r line
do
./scripts/bwa_aligner.sh $line
done<./lists/bwa_full_arg_list.txt
```

Once completed, have a look at the `bams` folder to make sure all our samples have been aligned.


#### Indexing BAMS 
Indexing aims to achieve fast retrieval of alignments overlapping a specified region without going through
the whole alignments. BAM must be sorted by the reference ID and then the leftmost coordinate before
indexing.
