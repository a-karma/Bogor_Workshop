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
| Symbol |Q-score|| Symbol |Q-score|
| ------ | ------ |------| ------ | ------ |
|!| 0||6|21|
|”|1||7|22|
|#|2||8|23|
|$|3||9|24|
|%|4||:|25|
|&|5||;|26|
|’|6||<|27|
|(|7||=|28|
|)|8||>|29|
|*|9||?|30|
|+|10||@|31|
|,|11||A|32|
|-|12||B|33|
|.|13||C|34|
|/|14||D|35|
|0|15||E|36|
|1|16||F|37|
|2|17||G|38|
|3|18||H|39|
|4|19||I|40|
|5|20||

> Question:  
> 
> Modify the command above to ...

Or even:
> Best practice: never use spaces in file names

note the empty line in the first quoted env (it seems to make it look nicer on github) 
