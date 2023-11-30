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

### Tutorial one - IQtree
IQTREE XXXXXX
read the docs XXXXXX

##### In this tutorial we will:
- use `screen`
- convert the files to phylip format 
- use IQtree to make an unrooted tree

Make session two directory 
conda env 
tutorial one 
screen 
vcf2phylip 
iqtree - without ASC 
iqtree - with ASC 
close screen 

### Tutorial two - exploring population structure with PCA
PCAs are used to XXXX

##### In this tutorial we will:
- get the files into the correct format
- run smartpca

Make sure you are still in the the correct conda environment for this session (the same as tutorial one). If not, reactivate with:

```sh
conda activate Day_3_b
```
Next make a new directory within session two for this tutorial. Something like `tutorial2_pca`
Enter this directory

##### 1. Convert files to the correct format
The first thing we need to do is convert our plink fileset into an eignstrat format which is used by the eigensoft set of programs - including smartpca. And to do this we use a program called `convertf`

Here we are going to generate a parameter file which contains the information that the program needs to correctly convert the files
> Make an empty text file and enter the text editor


Convertf can change between several different format, we have the .ped format and would like the EIGENSTRAT format as the output. 
> take a look at the format and additional options here - https://github.com/chrchang/eigensoft/blob/master/CONVERTF/README

So first you specific the location of your original plink ped and map file, define the outputformat that you require and then you want to generate the new eigenstrat files in your working directory.

```sh
genotypename:   /home/DATA/Day_3_b/babirusa_panel.ped
snpname:        /home/DATA/Day_3_b/babirusa_panel.map 
indivname:      /home/DATA/Day_3_b/babirusa_panel.ped
outputformat:    EIGENSTRAT
genotypeoutname: babirusa_panel.eigenstratgeno
snpoutname:      babirusa_panel.snp
indivoutname:    babirusa_panel.ind
```
> save with an informative name (e.g par.convertf_PEDtoEIGENSTRAT) and exit the editor

Now we can run convertf using 
```sh 
convertf -p [name_of_par_file]
```
Once completed, you should see the new files in your working directory

##### 2. Reassign the populations
Now check the first few lines of the `.ind` file. What do you see? The second column is sex (U = unknown) and the third is population. But because we dont want to make prior assumptions about the population membership of the individuals we want to rename this column so each individual is in a unique population. The easiest way is to copt the first column to the third. This can be done using `awk`

```sh
cat babirusa_panel.ind | awk 'BEGIN {OFS="\t"};{print $1,"U",$1}' > babirusa_panel.ind_new
```
We then need to remove the original `.ind` file and rename `.ind_new`
```sh 
rm babirusa_panel.ind
mv babirusa_panel.ind_new babirusa_panel.ind
```
##### 3. Run smartpca
Our data is now ready to run `smartpca`. Like `convertf` we need to make a parameter file to supply to the program to make the output files. 
> Make an empty text file and enter the text editor 

The same as the previous par file - we add the path of the input files followed by the output files
```sh 
genotypename:  babirusa_panel.eigenstratgeno
snpname:       babirusa_panel.snp
indivname:     babirusa_panel.ind
evecoutname:   babirusa_panel_PCA.evec
evaloutname:   babirusa_panel_PCA.eval
```
> save with an informative name (e.g par.smartpca_babirusa) and exit the editor

There are many other options we could add to the par file, which will depend on your data and the analysis you are conducting 
> Take a look at the documentation here: https://github.com/chrchang/eigensoft/blob/master/POPGEN/README

Now we run `smartpca` and save the output to a log file
```sh 
smartpca -p [name_of_your_par_file] > [name_of_your_logfile].log
```
This should give us two new output files:
- babirusa_panel_PCA.evec - this contains the XXXX
- babirusa_panel_PCA.eval - this contains the eigenvalues for each of the principal components (the importance of each axis)

In the next session we will copy these files to our local computer and visualise in RStudio. But now we will move on to the admixture analysis.

### Tutorial three - admixture analysis
Admixture XXX 

##### In this tutorial we will:
- run ADMIXTURE 

Again, make sure you are still in the the correct conda environment for this session (the same as tutorial one). If not, reactivate with:

```sh
conda activate Day_3_b
```
Next make a new directory within session two for this tutorial e.g `tutorial3_admixture`
Enter this directory
The program ADMIXTURE runs directly from a `.bed` file. We already have this in the right format.

### Lets quickly check up on the IQtree before having coffee...
Activate your screen session. To see which screens you have running, list them
```sh 
screen -ls
```
and then activate the one running IQtree
```sh 
screen -r [name_of_screen]
```
If it has not finished running leave it for longer while we have a break. Deattach the screen using `ctrl`+`a`+`d`
If it has finished - great :) 
When it has finished, check the files in your `tutorial1_tree` directory. You should have these output files: 
- file 1
- file 2
- file 3 ...


