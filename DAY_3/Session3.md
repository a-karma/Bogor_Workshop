![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day three - Session three - Visualisation 
In this session we are going to visualise all our results from the previous session. We will be working mainly on our local computers in RStudio and the browser.

We will:
- download our output files via `sftp`
- generate figures in RStudio for the PCA and ADMIXTURE results
- Visualise the phylogenetic tree in iTOL (https://itol.embl.de/) 

### 1. Downloading the output files 
On your local computer we will use `sftp` to connect to the remote server and download the files.

In your terminal make a project directory where you would like to store the output files on your local computer. 

Then run `sftp` - this stands for secure file transfer protocol and will open a connection between your local comput and the remote server. 
We need to tell it where the ssh key is, and which remote server to connect using your username and the IP address. Therefore enter: 
```sh
sftp -i ~/path/to/ssh/key [your_username]@138.246.238.65
```
You will now be in your home directory in the server, it should show this: 
```sh
sftp>
```

You are now connected to the server and can navigate through your directories like normal. 
If you add an `l` in front of your commands you should see you local computer file structure. Using the command `get` will download your files to your location on the local computer. 

You need to download your output files

- .evec and .eval for the PCA
- .Q for the admixture
- .treefile and the .iqtree files for the tree
  
using something like this:
```sh
ls
get tutorial_PCA/*ev* .
```

You will also need to download some files from the shared directory `/home/DATA/Day_3_b/`
- you need the population_file.txt
- and the rooted tree file `XXXX.treefile
When you have done this, type this to close the connection.

```sh
bye
```

### 2. Visualising and comparing the trees
You should have two trees now - the one the you make (unrooted) and the one from the shared directory (rooted).

The rooting of the tree requires an outgroup and to save time we did this for you. Briefly we used the same panel to call SNPs at the same location in the outgroup - which is the pygmy hog. This is then merged with the original fileset. The tree is reran in iqtree, and it is easy to specify this sample in the command with the option (`-o`).

Open your browser and navigate to: https://itol.embl.de/

Then open the .iqtree file in a text editor and copy the NEWICK format tree at the bottom of the file 
![NEWICK_tree](../IM/NEWICK_tree.png)

Paste it into the box on the site, and press upload. 

ADD THE INSTRUCTIONS FOR THE UNROOTING

ADD THE ROOTED TREE

WHAT ARE THE DIFFERENCES

### 3. Visualising the PCA analysis
Next we will look at the results of the PCA using RStudio.

Start by opening RStudio and open a new script
![RStudio_opening](../IM/RStudio_opening.png)

We will write our code into the script window (top left), where we can save everything we write. 
The code will run in the terminal (bottom left). You can also run commands here, but they will not be saved.
The environment on the top right is where all your data, objects and variables are stored. 
Plots are viewed on the bottom right. 

TYPE SOMETHING IN THE TERMINAL AND RUN - HOW DO YOU RUN ON A WINDOWS VS MAC

Ok so now, its good practice to clean your environment, set your working directory and load the correct libraries you need at the top of your script. Like this: 

```sh
rm(list=ls()) 
```
This will clean up any variables from a previous session which might conflict. 

We then need to set the working directory - this should be the path to the project folder you made, where you downloaded the output files to. 

```sh
setwd("~/PATH/TO/PROJECT/DIRECTORY/")
getwd() # this should show your working directory location in the terminal
```

Now you load the libraries you need for this section.
> These will have needed to be installed first with install.packages(), this only needs to be done once

```sh
library(tidyverse)
library(ggrepel)
```

Now we will load in the metadata file we need for PCA. This contains the sample names and the region of Sulawesi that the sample comes from. 

Load the metadata file and name the columns:
```sh
samplelist <- read_tsv("pop_file.txt", col_names = c("sample", "region"))
```
> `Hint` - if your file is in your working directory the you just type the name of the file. But you can use `tab` to navigate to where your file is. 

Check the top of your dataframe, you want to see the two columns
```sh
head(samplelist)
```

Next you want to make a new column which groups the samples from the north together.
First you make an empty column that is the length of the dataframe
```sh
samplelist$region2 <- rep(NA, length(samplelist$region))
```

Then you want to assign a new group based on the original grouping. So for the Southeast lets rename it to a group called "South" 
```sh
samplelist$region2[samplelist$region == "SE"] <- "South"
```
> Exercise
>
> Rename the Togean babirusa to a group called "Island" and then assign both the Northwest and Westcentral babirusa to a "North" group.

SAVE THE FILE FOR THE ADMIXTURE?

Next you will read in the datafiles, for the pca there are two files we need, `.evec` and `.eval`.
```sh
eval <- read.table("PATH/TO/FILE.eval")
```
> Exercise
>
> Can you read in the `.evec` file and check that they have be correctly loaded

Using the eigenvalues, we can calculate the percent contribution of each PC axis to the variation in the samples. This is for the first two axis.
```sh
evec.pc1 <- round(eval[1,1]/sum(eval)*100,digits=2)
evec.pc2 <- round(eval[2,1]/sum(eval)*100,digits=2)
```
> Extra exercise
>
> Can you calculate the percentages for the remaining 8 PCs?
> 

We will plot the PCA using ggplot and its easier to plot if everything is contained within the same dataframe. So we can add the metadata information to the eigenvectors. 

```sh
evec_merge <- as.data.frame(cbind(evec, sample = samplelist$sample, region = samplelist$region, region2 = samplelist$region2))
```
This command takes your eigenvalues, and adds the columns (`cbind()`) from the samplelist dataframe and names the columns. `as.data.frame()` ensures that the new dataframe is in the correct format.

Now we are ready for plotting using the function ggplot. ggplot is built in layers, with each layer being added to the basic plot using (`+`). 
The style of plot is determined by the `geom` functions. As we want a scatted plot we are using 

First lets make a basic scatter plot:
```sh
ggplot(data = evec_merge) + 
  geom_point(aes(x = V2, y = V3, colour = region), size = 4))
```
Lets unpack this
- the inital ggplot command determines the dataframe to be used
Then we need to specify the aesthetics (`aes`) within the geom function. These are the required values to build the plot correctly.
- the x-axis is specified as the second column which are the values for PC1
- the y-axis is specified as the third column which are the values for PC2
- the points will be coloured by the region factor
- size determines the size of the point

Lets make this an object called `pca_plot`

```sh
pca_plot <- ggplot(data = evec_merge) + 
  geom_point(aes(x = V2, y = V3, colour = region), size = 4))
```
And now to this we can add additional layers and make it look pretty
Lets add the percentage contribution to each axis
```sh
pca_plot <- pca_plot + xlab(paste0("PC1 (", evec.pc1, "%)")) + ylab(paste0("PC2 (", evec.pc2, "%)")) 
```
Finally lets make the plot look a little cleaner.
```sh
pca_plot <- pca_plot + theme_bw()
```
To see how it looks we need to just call the named object. 
ggplot has a huge amount of flexibility and there is alot of documentation online. You can also use `?ggplot` or `?geom_point` to read the help file.

![Final_PCA](../IM/Final_PCA.png)

Now when you look at the PCA, do you think that there is evidence of population structure? How many populations do you this these samples form?

You can save the plot with
```sh
ggsave(plot = PLOT_NAME, "PATH/TO/FILE.pdf")
```

### 4. Visualising the ADMIXTURE analysis
Based on these 






> Extra PCA exercises
>
> Can you add the sample names to the plot using the function `geom_text_repel()`
> Can you re-colour the points based the three clusters inside of four?
> Can you use scale_colour_manual() to choose your own colour palette?
> 




