![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
please do not modify these first two lines of the .md file

use the same syntax to add pictures:

placeholder name within square brackets and ../IM/file_name.png within parentheses

## Day two - Session two and three
In these sessions we are going to explore the population genetic structure and ancestry in the babirusa dataset. We assign the populations as XXXX
Explained in the lecture XXXX
In session two - we will run the anaysles on the virtual machine, building maximum likelihood tress in IQtree, a principal components analysis in smartpca and finally run ADMIXTURE with several values of K. 
In session three - we will export the data from the virtual machine and work in RStudio to visualise the results from the PCA and ADMIXTURE, and we will use the browser tool iTOL for tree visualisation.

### Tutorial two - exploring population structure with PCA
PCAs are used to

##### In this tutorial we will:
- get the files into the correct format
- run smartpca

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
