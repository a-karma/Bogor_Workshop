![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 3 - Population genetic analysis - Session 2
In this session we are going to learn the basics of population genomics analyses using the babirusa dataset. As explained, we will use these analyses to explore whether all the babirusa come from a panmictic population, or whether we can see evidence of population structure and whether this structure correlates with geography.

In this tutorial you will run the three analyses using the remote server:
- building maximum likelihood tress in IQtree
- principal components analysis in smartpca
- ADMIXTURE with different values of K

In session three - you will export the data from the virtual machine and work in RStudio to visualise the results from the PCA and ADMIXTURE, and you will use the browser tool iTOL for tree visualisation.

### 1. Building phylogenetic trees with IQtree

First create a new directory for this session in your home directory (`~/`). Recall the guidelines with went through on day 1.

Next activate the correct conda environment for the day
> `Hint` - this will be the same as the one you were using this morning

We will keep the output files for each analysis in their own directory, so make a directory for this analysis and navigate to it. This is your working directory.

### Convert to phylip format
The program you will use to make the tree is called `iqtree`. `iqtree` requires an alignment file as an input but this can take several formats. 
We will be using the `.phylip` alignment format. To generate this from the vcf file, you need to use a script called `vcf2phylip.py`. 

`iqtree` is program to reconstruct phylogenetic tree based on genetic data - it uses maximum likelihood (see here if you want more details about tree reconstruction methods: https://en.wikipedia.org/wiki/Computational_phylogenetics). 

The data for this exercise can be found here `/home/DATA/Day_3_b/`. 

First we need to to use `vcf2phylip.py` to convert our vcf (https://en.wikipedia.org/wiki/Variant_Call_Format) into a phylip (https://en.wikipedia.org/wiki/PHYLIP) format. This script can be foundin the scripts directory in the shared data folder. To do this, first create a shell variable we call `PANEL` which contains the path the panel.

```sh 
PANEL=/home/DATA/Day_3_b/babirusa_panel
```
Then we can run vcf2phylip by specifying the vcf file as the input (`-i`) and the output name you want in your current working directory (`--out-prefix`)

```sh
/home/DATA/Day_3_b/scripts/vcf2phylip.py -i $PANEL.vcf --output-prefix babirusa_panel
```
When it is finished you should see a your output file in your directory, and that with default settings the file has been named `[name_of_panel].min4.phy`

>  `Exercise one`
> 
>  Can you use the help file to find out what the addition of "min4" in the file name means? Now how would you run the command and change this option?

### Basics of using `screen` in linux
Now we have our input file in the correct format for `iqtree`. Because the command we are going to run can take a while to complete, we are going to use a program called `screen` to allow us to run commands in a different session in the background as we continue with the other analyses. Below are some of commands we will need for interacting with `screen`.

To activate new session with a specific name: 
```sh 
screen -S [name_of_session]
```
To detach a session you press `ctrl`+`a`+`d` at the same time 

To see what sessions are open use: 
```sh 
screen -ls
```

To reopen a session use: 
```sh 
screen -r [name_of_session]
```

To kill a session use: 
```sh 
screen -XS [name_of_session] quit
```

### Running iqtree
Now we have our input file and the basic understanding of screen we will run iqtree. 

Reopen the screen you made for running `iqtree`
You should see that you have to reactivate the conda environment when you enter a new session. 

IQtree is a very versatile program with many options yet very straightforward to use. 

To run it we are going to specify some parameters:
> 1. The phylip alignment (`-s`)
> 2. The sequence type (`-st`)
> 3. The substitution model (`-m`) - if you omit this option iqtree runs an automatic model test. But it can take a while so we will specify this. 
> 4. The output file prefix (`-pre`)

The theory behind substitutions model in phylogenetics is beyond the scope of this tutorial. We'd be happy to discuss this with you in the classroom if you are interested - you can also check this wikipedia page for more information: https://en.wikipedia.org/wiki/Substitution_model

So the command is: 
```sh 
iqtree -s babirusa_panel.min4.phy -st DNA -pre babirusa_panel_tree1 -m GTR+ASC
```

This should actually throw up an error and generate a new file for us. This is because `+ASC` model is specifying the ascertainment bias correction, which is appropriate for SNP data. ###LF COMMENT: WHY?### However there are invariable sites in the panel. The output file generated only contains variable sites. 

The reasoning behind the use of the `+ASC` model is beyond the scope of this practical - feel free to ask about this to an instructor if you want to know more.

> `Exercise two`
>
> How would you modify the original command to rerun iqtree using this new file?

If that is running, we will now let it in the background by deattaching the session (`ctrl`+`a`+`d`) and we will come back to it later.

### 2. Principal components analysis (PCA) in smartpca
PCA is a powerful statistical method widely used in population genomics to analyze and visualize patterns of genetic variation among individuals or populations. It helps researchers understand the underlying structure of genetic data and identify key trends that can be used to infer evolutionary relationships, migration patterns, and adaptation mechanisms.

In simple terms, PCA can be thought of as a technique that transforms a complex set of genetic data into a more manageable and interpretable form. It does this by identifying a set of new variables, called principal components (PCs), that capture the main axes of variation in the data. These PCs are essentially linear combinations of the original genetic markers, but they have the advantage of being uncorrelated with each other, making them easier to analyze and interpret.

Here we are going to use  `smartpca` to make a PCA from the babirusa data.  Make a new directory in your session two directory and navigate to it. 

### Convert files to the correct format
smartPCA does not work with vcf or plink files so we need to run a file conversion (this is somehow very common in bioinformatics). Luckily the author of the program provide a code to do the conversion from plink file to eigenformat. The first thing to do is convert our plink fileset into an eignstrat format which is used by the eigensoft set of programs - including smartpca. To do this we are going to use a program called `convertf`.

To run `convertf` you need to make a parameter file, or par file, which contains the information on where the files we want to convert are located. The program can change between several different formats, for example - the current files are in the .ped format and would like the EIGENSTRAT format as the output. 

> `Hint` - take a look at the format and additional options here - https://github.com/chrchang/eigensoft/blob/master/CONVERTF/README

Make an empty text file and enter the text editor (`nano`).

You need to specify the location of your original plink .ped and .map file, define the output format that you require and then specify where you want to generate the new eigenstrat files. For us this will be in your working directory.

```sh
genotypename:   /home/DATA/Day_3_b/babirusa_panel.ped
snpname:        /home/DATA/Day_3_b/babirusa_panel.map 
indivname:      /home/DATA/Day_3_b/babirusa_panel.ped
outputformat:    EIGENSTRAT
genotypeoutname: babirusa_panel.eigenstratgeno
snpoutname:      babirusa_panel.snp
indivoutname:    babirusa_panel.ind
```
> save with an informative name (e.g `par.convertf_PEDtoEIGENSTRAT`) and exit the editor

Now run convertf using 
```sh 
convertf -p [name_of_par_file]
```
Once completed, you should see the new files in your working directory

### Redefining the populations
Check the first few lines of the `.ind` file. What do you see? 
Column one should be the sample name, column two is sex (U = unknown) and column three is population. Currently the third colummn will be filled with `???`. We dont want to make prior assumptions about the population membership of the individuals, therefore we want to rename this column so each individual is in a unique population before we run the pca. The easiest way is to copy the first column to the third. This can be done using `awk`

```sh
cat babirusa_panel.ind | awk 'BEGIN {OFS="\t"};{print $1,"U",$1}' > babirusa_panel.ind_new
```

This command reads the .ind file with cat, this is then piped to awk which is printing the first field, followed by the "U" of the second field and then printing the first field again (i.e. the sample names). The output is written to a new .ind file. 

Then the original `.ind` file can be removed and rename the updated `.ind_new` file. Like this:
```sh 
rm babirusa_panel.ind
mv babirusa_panel.ind_new babirusa_panel.ind
```
These files can be used to run multiple types of analysis beyond PCA using the admixtools package (https://github.com/DReichLab/AdmixTools). We will stick to PCA for today though. 

### Run smartpca
Our data is now ready to run `smartpca`. Like `convertf` we make a parameter file to supply to the program to generate the output files.

Again, make an empty text file and enter the text editor.
The same as the previous par file - specify the path to the input files (in eigenstrat format generated by the convertf command), followed by the path we want for the two output files.

```sh 
genotypename:  babirusa_panel.eigenstratgeno
snpname:       babirusa_panel.snp
indivname:     babirusa_panel.ind
evecoutname:   babirusa_panel_PCA.evec
evaloutname:   babirusa_panel_PCA.eval
```
> save with an informative name (e.g `par.smartpca_babirusa`) and exit the editor

There are many other options to add to the par file, which will depend on your data and the analysis you are conducting 
> `Hint` take a look at the documentation here: https://github.com/chrchang/eigensoft/blob/master/POPGEN/README

Now run `smartpca` and save the output to a log file
```sh 
smartpca -p [name_of_your_par_file] > [name_of_your_logfile].log
```
This should give you two new output files, look at the contents of these files: 
- babirusa_panel_PCA.evec - this is the position of where each individual falls along the eigenvectors (columns 2:11) 
- babirusa_panel_PCA.eval - this contains the eigenvalues for each of the principal components (the importance of each axis)

> `Exercise three`
>
> Can you modify the par file and generate new outputs for 15 eigenvectors?

In the next session you will copy these files to your local computer and visualise in RStudio. But now lets move on to the admixture analysis.

### 3. ADMIXTURE analysis
The final analysis we will run is ADMIXTURE. ADMIXTURE is a software program for inferring individual ancestries and population structure from SNP data. It uses a maximum likelihood approach to estimate the proportions of each ancestral population in each individual. ADMIXTURE is a popular tool for population genomics research because it is relatively easy to use and can handle large datasets.

Make a new directory within your project folder for session two and navigate to it

The program `ADMIXTURE` runs directly from a `.bed` file. Great, we already have this in the right format. So the command is very simple.
```sh 
admixture [path_to_bed] [number_of_ks]
```
This program can also calculate the cross validation errors using the option (`-cv`).

> `Exercise four`
>
> Run ADMIXTURE for k = 2 with thr cross validation error calculation

Next you can use a loop to run through several values of k, calculate the cv errors and output the results to a log file.
First lets assign a variable for the path to the `.bed` file.

```sh
ADMIX=/home/DATA/Day_3_b/babirusa_panel.bed
```
Then make the for loop for k = 1-5
```sh 
for k in {1..5} 
do
    admixture $ADMIX --cv $k | tee babirusa_k_$k.log
done
```

When this has finished running, use `grep` on the log files to extract the cv values, which we can later visualise

```sh 
grep CV *.log > cv_errors.txt
```

### Lets quickly check on the IQtree before having coffee...

Reactivate your screen session. To see which screens you have running, list them first:
```sh 
screen -ls
```
and then activate the one running IQtree
```sh 
screen -r [name_of_screen]
```
If it has not finished running leave it for longer while we have a break. Deattach the screen using `ctrl`+`a`+`d`
If it has finished - great. You can check the files in the directory for this analysis. You should have these output files: 
- .iqtree - the report file
- .treefile - the ML tree in NEWICK format (this is the file we will need to visualise in the next session)
- .log - the log file
- .ckp.gz - this is a check point file, if a run is interrupted and you need to resume
- .mldist - maximum likelihood distance matrix
- .bionj - contains the BIONJ tree, related to the neighbour joining tree

If you see these files, you can close the screen session.

