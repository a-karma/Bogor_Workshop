![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 5 - Modelling gene flow 
In this session we will explore the `admixtool` suite and run some tests using f-statistics.

### Introduction

In day three you should have obtained these results from the PCA, iqtree and admixture analysis

![Workshop-logo](../IM/day5_recap.png)

Today we are going to test the support for the topology obtained from iqtree and whether RD7 is truely admixed

Before we start we should create the directory structure:
```sh
mkdir day5_fstat
cd day5_fstat
mkdir par;
ln -s /home/DATA/Day_5_s1/PANEL/ ~/day5_fstat/raw_data
```
In your `raw_data` folder you should see a plink file set consisting of two files, namely `chr10.ped` and `chr10.map`. In order to use them for our practical we need to convert them into eigenstrat format (like you did on day 3).

Like always - we will activate our conda env first
```sh
conda activate Day_5_s1
```

Let's start preparing the parameter file `~/day5_fstat/par/convertf.par` for conversion:
```sh
cd par
touch convertf.par
nano convertf.par
cd ..
```

Then inside the `convertf.par` file you need to include the following:
```sh
genotypename:  ~/day5_fstat/raw_data/chr10.ped
snpname:       ~/day5_fstat/raw_data/chr10.map 
indivname:     ~/day5_fstat/raw_data/chr10.ped
outputformat:    EIGENSTRAT
genotypeoutname: ~/day5_fstat/chr10.eigenstratgeno
snpoutname:      ~/day5_fstat/chr10.snp
indivoutname:    ~/day5_fstat/chr10.ind
```

You should now see the standard eigenstrat file set in your day5_fstat directory - three files `.eigenstratgeno`, `.snp`, `.ind`
Again, like day 3, we need to modify the `.ind` file  

```sh
cat chr10.ind | awk 'BEGIN {OFS="\t"};{print $1,"U",$1}' > chr10.ind_new
rm chr10.ind
mv chr10.ind_new chr10.ind
```

In order to calculate our f4 statistics, we need to prepare a configuration file for the program `qpDstat`. Which we will call `dstat.null.par` and should contain:
```sh
genotypeoutname: ~/day5_fstat/chr10.eigenstratgeno
snpoutname:      ~/day5_fstat/chr10.snp
indivoutname:    ~/day5_fstat/chr10.ind
popfilename:     ./par/null_test.txt
f4mode: YES
```

Finally we need to specify inside `./par/null_test.txt` file, which test we want to run. 
Therefore create this file inside the par directory - with the following contents
```sh
RD17 RD16 RD1 RD2
RD71 RD17 RD10 RD2
RD44 RD2 RD20 RD17
```
As you can see we are running three test (one per line) each involving 4 individuals (we are calculating f4s). The figure below shows the three tests we are specifying in the `./par/null_test.txt` file for the program to run.
The tests are:
1. two individuals from Togean (RD17, RD16) and two individuals from southeast(RD1, RD2) 
2. two individuals from Togean (RD71, RD17) and one north (RD10) and one southeast (RD2)
3. two individuals from the southeast (RD44, RD2) and one north (RD20) and one togean (RD17)
For each of the tests we assume that the topology obtained from iqtree is the correct one

![NULL](../IM/null-test.png)




