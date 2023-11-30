![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation
please do not modify these first two lines of the .md file

use the same syntax to add pictures:

placeholder name within square brackets and ../IM/file_name.png within parentheses

## Day two - Session two
In this session we are going to explore the population genetic structure and ancestry in the babirusa dataset.
We assign the populations as XXXX
Explained in the lecture XXXX
First, we will run all three analyses on the virtual machine, building maximum likelihood tress in IQtree, a principal components analysis in smartpca and finally run ADMIXTURE with several values of K. 
Then in session three - we will export the data from the virtual machine and work in RStudio to visualise the results from the PCA and ADMIXTURE, and we will use the browser tool iTOL for tree visualisation.

### Tutorial one - IQtree
IQTREE XXXXXX
read the docs XXXXXX

##### In this tutorial we will:
- convert the files to phylip format 
- use `screen` to run commands in the background
- use IQtree to make an unrooted tree

First in your home directory we need to make a new directory for this session. Call it something sensible and enter this directory. 

Next we will activate the correct conda environment for the day
> Hint - this will be the same as the one you were using this morning

Its good to keep all the output files for each analysis in their own directory. So make a directory for this tutorial (call it something like e.g. `tutorial1_tree`) and enter into it 

##### 1. Convert to phylip format
The program we will use to make the tree is called `iqtree`. This program needs an alignment as the input but this can take several format. 
We will be using the `.phylip` alignment format. To generate this, we need to use a script called `vcf2phylip.py`. You should be able to find this under `/home/DATA/Day_3_b/scripts/`. We want to use this script to convert our vcf panel into phylip format. So we can assign a variable the path to the panel, and supply this in the command with the `-i` option.

```sh 
PANEL=/home/DATA/Day_3_b/babirusa_panel
/home/DATA/Day_3_b/scripts/vcf2phylip.py -i $PANEL.vcf --output-prefix babirusa_panel
```
When it is finished you should see a new file in your directory, and that with default settings the file has been named `[name_of_panel].min4.phy`
This is because it applies a filter of a minimum number of samples per SNP and as default this is four individuals. You can change this using the option `-m`
> there are also some other options that can be changed, see https://github.com/edgardomortiz/vcf2phylip 

##### 2. Basics of using `screen` in linux
Now we have our input file for `iqtree`. Because the command we are going to run can take a while to complete, we are going to use a program called `screen` to allow it to run in the background as we continue with the other analyses. 
So to activate and name a new session we run: 
```sh 
screen -S [tree]
```
> To detach a session you press `ctrl`+`a`+`d` at the same time 
> To see what session are open you run: 
```sh 
screen -ls
```
> To reopen a session run: 
```sh 
screen -r [name_of_session]
```
> To kill a session run: 
```sh 
screen -XS [name_of_session] quit
```

##### 3. Running iqtree
Now we have the basics. Reopen the screen you made for running `iqtree`
You should see that you have to reactivate the conda environment as well, when you enter a new session. 

IQtree is a very versatile program with many options but is very simple to use. 
To run it we are going to specify some parameters:
> 1. The phylip alignment (`-s`)
> 2. The sequence type (`-st`)
> 3. The substitution model (`-m`) - if you omit this option iqtree runs an automatic model test. But it can take a while so we will specify this. 
> 4. The output file prefix (`-pre`)

So the command is: 
```sh 
iqtree -s babirusa_panel.min4.phy -st DNA -pre babirusa_panel_tree1 -m GTR+ASC
```

This should actually throw up an error and generate a new file for us. This is because `+ASC` model is specifiying the ascertainment bias correction, which is appropriate for SNP data. However there are invarible sites in the panel. The output file generated only contains variable sites. 

> How would you modify the original command to rerun iqtree using this new file? "babirusa_panel_tree1.varsites.phy"

If that is running well. You can let it run in the background by deattaching the session (`ctrl`+`a`+`d`). We will come back to it later.

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

Again, make sure you are still in the the correct conda environment for this session (the same as tutorial one and two). If not, reactivate it.

Next make a new directory within session two for this tutorial e.g `tutorial3_admixture`
Enter this directory
The program ADMIXTURE runs directly from a `.bed` file. We already have this in the right format.

First lets assign a variable for the path to the `.bed` file 
```sh
ADMIX=/home/DATA/Day_3_b/babirusa_panel.bed
```

Then the program is very simple its just
```sh 
admixture [path_to_bed] [number_of_ks]
```
Can you run this for k=2? 

However we can also calculate the cross validation errors (cv). These are XXXXXX
To do this we will set up a loop to run through several values of k while calculating the cv errors and outputting it to a log file. 
```sh 
for k in {1..5} 
do
    admixture $ADMIX --cv $k | tee babirusa_k_$k.log
done
```

When this has finished running, we can use `grep` on the log files to extract the cv values, which we can later visualise

```sh 
grep CV *.log > cv_errors.txt
```

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


