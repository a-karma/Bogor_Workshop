![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

##  

Fastq Format
Genomic information from high-throughput sequencing are stored in text-based files called
FASTQ. These files comprise a series of entries containing not only the sequence but also
its quality score. Each entry consists of four lines as shown in the box below.

```sh
@HiSeq 4000:1:FCX:4:15:66:165 2:N:0:7
ATTTAGTACCATGACATGACACATACTACAATTGACGACATCAATCA
+
IGHFDEC@;;?=>B=?<;A:?@>9<>9756867544312*+*),&)’
```

The first line starts with @ and contains a series of characters which uniquely identify
the read. When dealing with Illumina data, the identifier includes 7 fields encoding
information about the sequencing process plus 4 fields about the read itself (see tab2.1).
Fileds are separated by a : and the two groups are separated by a space.
Table 2.1: Fastq format – Fields in the Sequence ID line
Field
Instrument used for sequencing
Run Number on that instrument
flow cell ID
lane number
tile number
X coordinate of cluster
Y coordinate of cluster
read number
is filtered
control number
sample number
Frantz-Lab, QMUL, 2020
Example
HiSeq 4000
1
FCX
4
15
66
165
2 (2 nd read of the pair)
N (Y did not pass, N otherwise)
always zero on HiSeq X
7 (the 7 th sample in the list)

### Example of section title 
##### Example of sub-section title 
- this is how you define bullet point
- in case you want to make more explicit 
- the series of steps required to accomplish a task 

See example below on how to format commands that the participants will have to run

```sh
conda activate Day_1
plink --bfile file_name --recode
```

you can instead use `this syntax` to highlight an in-line command, software name or something you think it's important

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
