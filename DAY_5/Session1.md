![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 5 - Modelling gene flow 
In this session we will explore the `admixtool` suite and run some tests using f-statistics.

### Introduction

In day three you should have obtained these results from the PCA, iqtree and admixture analysis


![Workshop-logo](../IM/day5_recap.png)




Before we start we should create the directory stracture:
```sh
mkdir day5_fstat
cd day5_fstat
mkdir par;
ln -s /home/DATA/Day_5_s1/PANEL/ raw_data
```
In your `raw_data` folder you should see a plink file set consisting of two files, namely `chr10.ped` and `chr10.map`. In order to use them for our practical we need to convert them into eigenstrat format (like you did on day 3).

Let's start preparing the parameter file `~/day5_fstat/par/convertf.par` for conversion:
```sh
content here
```
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
