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

The `GTF` and the `BED` format are both TAB-separated files used to store genomic regions as coordinates and their associated annotations. 
Although in principles we could manually edit these files using standard text editors this becomes very unpractical when dealing with very large files.
Alternatively, we could combine command line tools (like `sed`, `awk`, or `grep`) to make this process more efficient but performing complex tasks using only these tools is not straightforward and often require a deep understanding of programming. Luckly for us, bioinformaticians have created various software specifically designed to manipulate these file formats. 

Let's have a look at `bedtools` a powerful toolset for genome arithmetic. 

In your terminal, please type: 

```sh
bedtools --help
```

Looks like the program is not installed :( 

On servers you often do not have permissions to install softwares directly. This is because 
You can read more about conda [here](https://docs.conda.io/en/latest/)



```sh
conda activate Day_1
plink --bfile file_name --recode
```

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
