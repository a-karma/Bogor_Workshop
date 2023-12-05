![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 2 - NGS data pre-procesing - Session 1

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
Fields are separated by a : and the two groups are separated by a space.

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



you can also use this env for exercises and tips:
> Exercise 1 
> 
> Modify the command above to ...

Or even:
> Best practice: never use spaces in file names

note the empty line in the first quoted env (it seems to make it look nicer on github) 
