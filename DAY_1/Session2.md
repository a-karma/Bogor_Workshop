![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
## Day 1 - Session 2: Basic concepts of command line programming

### 1. Shell Scripting
In Session 1 we have seen how to navigate a Unix-like file system and how to manipulate text files. 
In this section we will revise what we have learnt about variables in Bash, we will then introduce the concept of positional arguments, and we will write our first shell script.

Before we start, let's run some preliminary commads to create the directory structure for this session:

```sh
cd; mkdir session2; cd session2; mkdir raw_data; mkdir scripts; mkdir results 
```
Note the use of `;`, which allows to run multiple commands in short succession.  

> Exercise 1
>
> Create a symbolic link between the `/home/DATA/Day_1` folder and your newly created `raw_data` directory 

If you now move into your raw_data directory and run `ls Day_1` you should see two files having the `.txt` extension, namely:
```sh
instructors_list.txt
participants_list.txt
```

Let's have a look at `instructors_list.txt` first, you can print the content on screen using `cat`.

Unfortunately, the fields are not well defined because each word is separated by a space. let's try to fix these formatting issues.

First of all we need to separate the last two fields (affiliation and status) from the instructors' names.
We can do this in `awk` and make use of variable `NF` (number of fields variable) which is set to the total number of fields in the input record:
```sh
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/Day_1/instructors_list.txt
```
Now let's redirect the output to a file:
```sh
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/Day_1/instructors_list.txt > ~/session2/raw_data/aff_status.txt
```
Tabs allow for consistent indentation across different file viewers and analysis tools. Each tab represents a single logical level of indentation, making the alignment structure more apparent and easier to interpret. Tabs are also more compatible with automation tools and scripts used for sequence analysis and manipulation. Many software programs and scripts expect tab-delimited alignment files, and using tabs avoids potential compatibility issues or the need for manual formatting adjustments.

Let's now deal with the names. Given that we don't know how many words each names consit of, we should start by printing all but the last two fields of the original input file:

```sh
awk 'NF-=2 {print $0}'  ~/session2/raw_data/Day_1/instructors_list.txt
```

Then we can pipe this into `sed` and replace all white spaces with the character `_`:

```
awk 'NF-=2 {print $0}'  ~/session2/raw_data/Day_1/instructors_list.txt | sed -e 's/ /_/g'
```
Note the use of the `g` at the end of the substitution command. The use of `g` means that sed will now change all matching space into a _.

Finally we redirect the output to a file:

```
awk 'NF-=2 {print $0}'  ~/session2/raw_data/Day_1/instructors_list.txt | sed -e 's/ /_/g' > ~/session2/raw_data/names.txt
```

Now that we have created these two intermediate files we can stitch them together to reconstruct the initial information correctly formatted:
```sh
cd ~/session2/raw_data/
paste names.txt aff_status.txt > corrected_instructors_list.tsv
```
The `paste` command is a versatile tool that combines text from multiple files. It's particularly useful for merging data, and comparing files.

The series of commands presented above acts as a single unit to accomplish the required task. Therefore we can transform them into a script that we can re-use.
```sh
cd ..
touch ./scripts/formatting.sh
```
Note the `cd ..` command which moves the user up one directory (i.e. from `raw_data` to `session2` in this case).

Now we need to transform into an executable file. to do so run:
```sh
chmod 770 ./scripts/formatting.sh
```
In unix-like systms, `chmod` is the command and used to change the access permissions. In this case the owner of the file (i.e. you) should now be able to read (4), write (2), and execute (1) this file, hence the first 7 which is the sum of the three permission granted. The same is true for other users in the group (the second 7) while external user do not have any permission (0). 

Now let's use nano to edit our script and add this code-block to its content:

```sh
#!/usr/bin/bash
awk '{print $(NF-1),"\t",$NF}' ~/session2/raw_data/Day_1/instructors_list.txt > aff_status.txt
awk 'NF-=2 {print $0}' ~/session2/raw_data/Day_1/instructors_list.txt | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt >  ~/session2/results/corrected_instructors_list.tsv
```

This version of the script formatting.sh is not very useful because it can work only on the `instructor_list.txt` input file.
If we want to re-use it to format a different input file we would have to open it and edit the file name every time which is not convenient.
Let's modify it to allow for more flexibility in the usage by transforming input and output into variables:

```sh
#!/usr/bin/bash

INPUT_FILE=~/session2/raw_data/Day_1/instructors_list.txt
OUTPUT_FILE=~/session2/results/corrected_instructors_list.txt

awk '{print $(NF-1),"\t",$NF}' $INPUT_FILE > aff_status.txt
awk 'NF-=2 {print $0}' $INPUT_FILE | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt > $OUTPUT_FILE
rm names.txt
rm aff_status.txt
```
Note that we have added two lines to delete the intermediate files (names.txt and aff_status.txt) that we don't need anymore.

This version looks slightly better but the input and output are still hard-coded inside the script. 
Ideally, we would like to supply the input and output at the runtime (meaning when we execute the script). To do so we can make use of positional arguments.
The indexing of the arguments starts at one, and the first argument can be accessed inside the script using $1. Similarly, the second argument can be accessed using $2, and so on.
Thus our final version of `formatting.sh` should be:
```sh
#!/usr/bin/bash

INPUT_FILE=$1
OUTPUT_FILE=$2

awk '{print $(NF-1),"\t",$NF}' $INPUT_FILE > aff_status.txt
awk 'NF-=2 {print $0}' $INPUT_FILE | sed -e 's/ /_/g' > names.txt
paste names.txt aff_status.txt > $OUTPUT_FILE
rm names.txt
rm aff_status.txt
```
now we can execute the script from the `session 2` directory and provide the correct input and output at the call:

```sh
./scripts/formatting.sh ./raw_data/Day_1/instructors_list.txt ./results/corrected_instructors_list.txt
```

> Exercise 1
>
> Use this script to correctly format the `participants_list.txt`. Make sure you input the right path to this file and use a sensible name for the output file (second argument)
 

### 2. Working with bioinformatic softwares using conda
In Session 1 we have seen three examples of text file that are commonly used in bioinformatics:

- The `fasta` format to store DNA sequence information
- The `GTF` (Gene Transfer Format) developed for the Ensembl genome browser
- The `BED` (Browser Extensible Data) format developed at UCSC for the Genome Browser tool

The `GTF` and the `BED` format are both TAB-separated files used to store genomic regions as coordinates along with their associated annotations (see the previous session for more information on `GTF`). While `GTF` is used to store mostly information about genic structure accross the genome (i.e. the G stands for gene), the `BED` (short for Browser Extensible Data) stores information about any genomic features. It is simpler than `GTF` and widely used in bioinformatics and consist of at least 3 fields:

1. chromsome ID
2. start of feature (0 based coordinate; see link to UCSC below)
3. end of the feature 

Other fields (4 etc.) are free form (as long as they are tab separated) - these can be used to store any information you may want, e.g. the name of the interval. For more information about this widely use format refer to the UCSC website here: https://genome.ucsc.edu/FAQ/FAQformat.html#format1

Although in principles we could manually edit these files using standard text editors this becomes very unpractical when dealing with very large files.
A more efficient option would be combinig command line tools (like `sed`, `awk`, or `grep`) but performing complex tasks using only these tools is not straightforward and often require a good basis understanding of unix. Luckly for us, bioinformaticians have created various software specifically designed to manipulate these file formats. 

Let's have a look at `bedtools` a powerful toolset for genome arithmetic that focuses on `BED` format.

In your terminal, please type: 

```sh
bedtools --help
```

Looks like the program is not installed :( 

To protect the integrity of the file system on a server, normal users do not have permissions to directly install softwares. Moreover, almost any bioinformatic tool will rely on specific libraries or package versions that might create conflicts or even impeed the functionality of other programs. To circumvent these issues, we need to make sure that the software we need are installed in a "confined space" (a.k.a environment) containing all the necessary dependencies which becomes accessible only when we need it. This can be easily implemented using `conda` which is a package, dependency, and environment manager for any programming language. 

You can read more about conda [here](https://docs.conda.io/en/latest/).

We have already installed conda on our server and we have created a different environment for each day of the workshop. In order to have access to conda please run:

```sh
source /home/anaconda3/bin/activate
conda init
```
You need to run this command only once and you should see a change in your prompt: the word `(base)` appears on the left.
This is signalling that you are now in the conda base environment which represent the default space. 
To activate the environment for this session, simply run:

```sh
conda activate Day_1
```

> Question: Have a look at your prompt again, what do you see? 

After activating an environment all software installed in it become immediately accessible. Let's check whether we can use bedtools now:

```sh
bedtools --help
```

Hurray! Now that bedtools is accessible let's see what we can do with it.

The `snp_ch30.bed` file in the folder `/home/DATA/Day_1/` is an example of a "customized" bed format. It contains the three mandatory fields (chromosome, start, end) plus an unusual 4th field. In that column I have stored the genotype of 4 individuals at that position. If the 4 th column has a lowercase `m`, that particular site is monomorphic.
You can have a look at it using `less` or inspect just three lines with a combination of head and tail, see for example what you get by running:

```
head -65 ~/session2/raw_data/Day_1/snp_ch30.bed | tail -3
```

Suppose you are interested in analysing neutral evolving sites, therefore, you may want to remove from the analysis all sites that are likely to be under selective pressures. 
As a first approximation, we could take a conservative approach and start to analyse polymorphic sites (SNPs) that are not in coding region (i.e. do not code for protein and so less likely to be under selective constraints; please speak to your instructor if you do not understand this concept). Performing this task manually is obviously tedious and very time consuming but it's super fast using a software like bedtools:

```sh
cd ~/session2/
bedtools intersect -a ./raw_data/Day_1/snp_ch30.bed -b ./raw_data/Day_1/genes_chr30.gtf -v > ./results/snp_filtered.bed
```
The intersect command reports overlapping regions between two BED/GFF/GTF files by comparing the coordinates of the genomic feature listed in them.
The `-v` flag tells `intersect` to report all lines in file A (specified using the `-a` flag) that DO NOT overlap with the genomic intervals listed in file B (-b flag).

> Exercise 2
>
> How many SNPs we have excluded?
>
> hint: remember that each SNP information is recorded on a single line in the bed file format

There is a lot more that you can do with genome arithmetic. Let’s picture a more complex scenario: suppose you are interested in studying the promoter regions of various genes on this chromosome. You have a fasta file with the entire chromosome sequence (see `ptw_ch30.fa` file) and you would like to examine 10 kb upstream the starting codon of each gene. 

How can you get that information?

First of all you need to get the coordinates of the starting codons. This is very easy using grep:
```sh
grep 'start_codon' ./raw_data/Day_1/genes_chr30.gtf > ./raw_data/CDS_start.gtf
```

Now you need to modify this file using the function `flank` implemented in bedtools which will create flanking intervals for each region in a BED/GFF/VCF file.
This function requires also a genome file defining the length of each chromosome, so let’s create this file first.
```sh
echo -e "chr30\t40214260" > ./raw_data/ch30_length.bed
```
The above command `echo` would normally output on screen whatever string you type after it. The value `40214260` is the size of chromosome 30 in the dog genome (which we used for this example).
The `-e` option tells the software to enable the interpretation of backslash escapes and `\t` stands for TAB.
Now we can run:

```sh
bedtools flank -i ./raw_data/CDS_start.gtf -g ./raw_data/ch30_length.bed -l 10000 -r 0 > ./raw_data/ch30_promoters.gtf
```
where the `-l` and `-r` flags stand for left (or upstream) and right (or downstream) and the number following each of these flags represents the length of the flanking region.

Finally we can use the `getfasta` function in bedtools to extract the actual sequences of the promoter regions:

```sh
bedtools getfasta -fi ./raw_data/Day_1/ptw_ch30.fa -bed ./raw_data/ch30_promoters.gtf -fo ./results/ptw_prom_sequences.fa
```

where the `-fi` flag stands for file input, the `-bed` indicates the coordinate file while the `-fo` option stands for file output. 

You can examine the first output line using: 
```sh
cat ./results/ptw_prom_sequences.fa | head -1
```

> Exercise 3
>
> Combine the intersect and flank functions in order to filter the `snp_ch30.bed` file
> by excluding CDS and all regions that are 5 kb from the start and the stop codon of each CDS.

### 3. Transferring files
We conclude today's tutorials with some basic commands to transfer files between you local machime and our server.
You can download or upload files using the Secure File Transfer Protocol (SFTP) which essentially creates a secure connection between your laptop and the cloud server.

##### Windows Users
In the guide we've sent you, you've seen how to install putty and how to connect to our cloud server using this program. 
You are going to use a different component of putty to open an sftp connection. In the search bar please type 'PSFTP'.
 - open the PSFTP app like in figure
 - In the terminal window type `open Bogor_Workshop` 
 - You should be now asked to insert the log-in name, please insert your user_ID
 - Finally you should be asked for your key-passphrase, please type the same one you have use for ssh connection and press enter

![psftp](../IM/psftp.png)

##### Linux & macOS user
Please open a terminal and then run:

```sh
sftp -i ./.ssh/workshop_key user_ID@138.246.238.65
```
where user_ID needs to be replaced by the credentials we sent you.

##### Commonly used sftp commands
After establyshing an sftp connection you can use `lcd` to change directory on your local machine while you can use the standard `cd` command
to navigate the cloud server file system. 

In order to upload a file on the server simply type:
```sh
put file_name
```
If instead you want to download a file from the server use:
```
get file_name
```
To close the sftp connection, simply type `bye`

If you finished early and you want to keep playing with bedtools have a look at this tutorial here: https://sandbox.bio/tutorials?id=bedtools-intro

