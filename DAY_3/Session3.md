![Workshop-logo](../IM/LOGO_new.png)
# Applications of Genomics in Wildlife Conservation

## Day 3 - Population genetic analysis - Session 3 
The aim of this session is to visualise the results of the analyses we just completed in the previous session. We will be working mainly on our local computers in RStudio and a web browser.

Before we can do that we first need to download our output files via `sftp`. We will then:
1. visualise the phylogenetic tree in iTOL (https://itol.embl.de/)
2. generate figures in RStudio for the PCA and ADMIXTURE results

### 1. Downloading the output files 
On your local computer please use either the `sftp` command or the `PSFTP` app to connect to the remote server and download files.

Please create a project directory where you would like to store the output files on your local computer. Like you have done on the previous days but now on your local machine. Make a directory for the output of each analysis and one directory for storing the figures

For example:
```sh
mkdir day3; cd day3; mkdir tree_output; mkdir pca_output; mkdir admixture_output; mkdir figures
```

##### macOS, Linux, and WSL users
Open a terminal and please enter: 
```sh
sftp -i ~/path/to/ssh/key [your_username]@138.246.238.65
```
You will now be in your home directory in the server. 

##### Windows users
Open the PSFTP app then type:
```sh
open name_of_the_putty_saved_session
```
Then enter your user_ID then the key passphrase.

You will now be in your home directory in the server.

Either case, you should now see a terminal like the one in the screenshot below. Here you can see: 
- the command `ls` should list the files in your home directory
- If you add an `l` infront, so `lls` you should see the file structure for your local computer
- `pwd` shows the working directory in the server
- while `lpwd` shows the working directory on your local computer
 
![SFTP_connection](../IM/SFTP_connection.png)

We need to download the output files from the previous analyses:

- From the PCA = `.evec` and `.eval`
- From the ADMIXTURE = all the `.Q` files and `cv_errors.txt`
- For the trees = `.treefile` and the `.iqtree`

Using the command `get` will download your files to your location on the local computer. 
The structure is like this: 
 ```sh
get [PATH/TO/REMOTE/FILE] [LOCATION/TO/PUT/FILE/IN/LOCAL]
```

> Hints - You can use `get -r ` to get the whole directory.
> The -r flag, also known as --recursive, is a powerful option used in various shell commands to operate on files and directories in a nested manner. It instructs the command to process not only the specified target but also all of its contents, including subdirectories and their contents, recursively.
> 
> If you are in the desired location on your local you can just use a `.` instead of `[LOCATION/TO/PUT/FILE/IN/LOCAL]`

You will also need to download some files from the shared directory `/home/DATA/Day_3_b/`
- you need the `pop_file.txt`
- and the rooted tree file `babirusa_rooted.treefile`
  
When you have downloaded all the files, type this to close the connection.
```sh
bye
```

Now move the output files to the correct directory for each analysis.

### 2. Visualising and comparing the trees
You should have two trees now - the one the you built (unrooted) and the one from the shared directory (rooted).

The tree are in a newick format, a text-based format used to represent phylogenetic trees. It is widely used in bioinformatics and evolutionary biology to store and share tree data. For more information check out this wikipedia: https://en.wikipedia.org/wiki/Newick_format

The rooting of the tree requires an outgroup and to save time we did this for you. Briefly we added the genotype of an outgroup (pygmy hog) to the SNP panel before computing the tree.

Open your browser and navigate to: https://itol.embl.de/

Scroll down a little and click on the upload tree button

#### Unrooted tree
First we will look at the unrooted tree you generated. So under the tree file box, navigate to your unrooted `.treefile` and click upload. Alternatively you can copy the newick format tree at the bottom of the `.iqtree` file into the text box.

> `Exercise one`
>
> Use the control panel on the right to make your unrooted tree look like this example below

![iTOL_unrooted_renamed](../IM/iTOL_unrooted_renamed.png)

#### Rooted tree
Now lets do the same for the rooted tree, in a different browser window so you can compare. The first thing to do is specify the root. Then we will colour the samples by region to identify if the populations form monophyletic clades. 

> `Exercise two`
>
> By hovering over the samples, root the tree by the correct sample and colour by region using the `pop_file.txt`

It is possible to add colours for regions, change and format lines etc, but for speed we will just colour the groups based on the population metadata file ('pop_file.txt') so we can compare the clades we see on our tree and see if there is any geographic signal. 

> `Exercise three`
>
> Using the population metadata, and the `coloured ranges` options in iTOL can you colour the samples by region. The image below shows the Togean island population. Can you do the same for the individuals from the SE, WC and NW. But you might need to consider how to colour the samples from the NW and WC? 

![iTOL_rooted_TO](../IM/iTOL_rooted_colours_TO.png)

#### Questions:
- when you consider the region of origin for each individual, do you think there is evidence of geographic clustering?
- what is the difference between the tree with and without the outgroup?
- which population of babirusa are most distinct?

### 3. Visualising the PCA analysis
Next we will look at the results of the PCA using RStudio.

Start by opening RStudio and open a new script (pink arrow in the screenshot).
![RStudio_opening](../IM/RStudio_opening.png)

- we will write our code into the source window (top left), where we can edit and save everything we write. 
- the code will run in the terminal (bottom left). You can also run commands here, but they will not be saved.
- the environment on the top right is where all your data, objects and variables are stored. 
- plots etc are viewed on the bottom right. 

When we write something in the script, it is not run until you source it. 
First highlight the code you want to run then:
- mac - press `cmd` + `entr`
- windows - press `ctrl` + `entr`

Ok so now, its good practice to clean your environment, set your working directory and load the correct libraries you need at the top of your script. Like this: 

```sh
rm(list=ls()) 
```
This will clean up any variables from a previous session which might conflict. 

We then need to set the working directory - this should be the path to the project folder you made, where you downloaded the output files to. 

```sh
setwd("~/PATH/TO/PROJECT/DIRECTORY/") # this would be the `day3` project directory
getwd() # this should show your working directory location in the terminal
```
You can see above that anything your write with a `#` will not be run - it is a comment. In programming, a commented line is a line of code that is intentionally excluded from execution. Comments are typically used to explain the purpose of code, document functionality, or provide additional information for programmers who may need to understand or modify the code in the future.

Now you load the libraries you need for this section.
> Hint - These will have needed to be installed first with `install.packages()`, and this only needs to be done once

```sh
library(tidyverse)
library(ggrepel)
library(ggpubr)
```
In R programming, libraries are collections of functions and data sets that provide additional functionality and capabilities for data analysis, visualization, and other tasks. Libraries are essential for extending the abilities of the base R environment and making it a more powerful and versatile tool for various applications.

Now we will load in the metadata file we need for plotting the PCA, which contains the sample names and the region of Sulawesi that the sample comes from. 

Load the metadata file and name the columns:
```sh
samplelist <- read_tsv("pop_file.txt", col_names = c("sample", "region"))
```

In R, the `read_tsv()` function is used to import tab-separated values (TSV) files into R data frames. TSV files are a common format for storing tabular data, with each line representing a record and each column representing a variable. The columns are separated by tab characters, hence the name "tab-separated values".

> `Hint` - if your file is in your working directory the you just type the name of the file. But you can use `tab` to navigate to where your file is.

Check the top of your dataframe, you want to see the two columns
```sh
head(samplelist)
```

Next you want to make a new column which groups the samples from the north together.
First you make an empty column that is the length of the dataframe. In R `$` indicates the column name you want to specify. RStudio is helpful and when you type `$` it will show you a list of the potential columns in that dataframe.
```sh
samplelist$region2 <- rep(NA, length(samplelist$region))
```

Then you want to assign a new group based on the original grouping. So for the Southeast lets rename it to a group called "South" 
```sh
samplelist$region2[samplelist$region == "SE"] <- "South"
```
> `Exercise three`
>
> Rename the Togean babirusa to a group called "Island" and then assign both the Northwest and Westcentral babirusa to a "North" group.

Next you will read in the datafiles, for the pca there are two files we need, `.evec` and `.eval`.
```sh
eval <- read.table("pca_output/[name_of_file].eval")
head(eval)
```
> `Exercise four`
>
> Can you read in the `.evec` file and check that it has been correctly loaded - there should be 12 columns

Using the eigenvalues, we can calculate the percent contribution of each PC axis to the variation in the samples. This is for the first two axis.
```sh
evec.pc1 <- round(eval[1,1]/sum(eval)*100,digits=2)
evec.pc2 <- round(eval[2,1]/sum(eval)*100,digits=2)
```

We will plot the PCA using ggplot and first we will add the metadata information to the eigenvectors dataframe.  
```sh
evec_merge <- as.data.frame(cbind(evec[,-12],
                                  sample = samplelist$sample,
                                  region = samplelist$region,
                                  region2 = samplelist$region2))
```
- `cbind()` function binds the columns together
- we take all but the last column of the eigenvector df (this was an additional sample name column)
- plus the metadata columns we want - the format is `[name_of_column] = [location_of_column]`
- `as.data.frame()` function ensures that it stays in a data frame format

> `Question`
>
> Why do we do `evec[,-12]` first? What is it doing? `Hint` - check the contents of the original `evec` file

Now we are ready for plotting using the function ggplot. ggplot is a very comprehensive package with many options. We will just try and keep it simple but if you have experience with it, go ahead and add options ontop. 
The plots are built in layers, with each layer being added to the basic plot using (`+`). 
The style of plot is determined by the `geom` functions - for a scatter plot we want to use `geom_point()`

First lets make a basic scatter plot:
```sh
ggplot(data = evec_merge) + 
  geom_point(aes(x = V2, y = V3, colour = region), size = 4)
```

Lets unpack this
- the inital ggplot command determines the dataframe (`data = evec_merge`) to be used
- the `geom_point()` specifies the point or scatter plot we want
  
Then we specify the aesthetics (`aes`) within the geom function. These are the necessary values to build the plot correctly.
- the x-axis is specified as the second column which are the values for PC1
- the y-axis is specified as the third column which are the values for PC2
- the points will be coloured by the region factor - this is important so we can see if there are any distinctive clusters
- size determines the size of the point

ggplot2 is a powerful data visualization library for the R programming language see here for more detials on ggplot2: https://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html

We will do the same thing again but this time specifying an object name. Lets make this an object called `pca_plot`.
```sh
pca_plot <- ggplot(data = evec_merge) + 
  geom_point(aes(x = V2, y = V3, colour = region), size = 4)

pca_plot
```
Then you can see the plot by calling the object `pca_plot` on its own

Now this your basic PCA plot. You could stop here - but there are many other options that we could add, and just adding a couple can make the figure much clearer for the reader.

First lets add the percentage contribution to each axis, this gives an indication of how important each axis is to the genetic structure.
```sh
pca_plot <- pca_plot + xlab(paste0("PC1 (", evec.pc1, "%)")) + ylab(paste0("PC2 (", evec.pc2, "%)")) 
```

Next lets make the plot look a little cleaner. ggplot has lots of built in `theme` options.
```sh
pca_plot <- pca_plot + theme_bw()
```
> `Hint` - start typing `theme_` and then use the tab to see the different options

The theme_bw() function is part of the ggplot2 package in R and is used to create a black and white theme for your plots. 

You are saving the layer ontop of the initial object, so be warned if you want to make changes or clear something you will have to run the code from the beginning.

Finally - we will add the sample names. This is useful so that we can identify any outliers in the data or interesting individuals. Because this needs a new geom, we will give this to a new object so we do not confuse R.

```sh
pca_plot_names <- pca_plot + geom_text_repel(aes(x = V2, y = V3, colour = region, label=sample))

pca_plot_names
```

<img src="../IM/Final_PCA.png" width="80%" height="80%">

> `Hint` - use `?ggplot` or `?geom_point` to read the help files and there is alot of help online

You save the plot using `ggsave`. Save the PCA as a `.png` in your figures directory.
```sh
ggsave(plot = PLOT_NAME, "figures/name_of_plot.png") # you can make a .pdf or other formats by changing the extension here
```

#### Questions:
- when you look at the PCA, what do you see and what do you think this means for the babirusa population?
- do you think that there is evidence of population structure and how many populations do you this these samples from?

### 4. Visualising ADMIXTURE analysis
Great, now lets move on to plotting the admixture results. You can keep going in the script, but remember to save it every so often. 
We want to use the same metadata file as for the PCA, so lets not clear our environment. 

First we need to read in the `.Q` files we downloaded from the server. 
```sh
k2 <- read.table("PATH/TO/FILE.Q")
head(k2) # you should see two columns "V1" and "V2" 
```
The columns V1 and V2 correspond to the ancestry proportions, and the number of columns will change with the values of K
> `Exercise five`
>
> Can you read in the additional Q files for the each of the clusters you ran?

Again, the first thing we need to do is get the dataframe in the correct format. We can use the `cbind()` function from before, or because we also need to get the data into longform format, we can achieve this using the `tidyverse` library (of which ggplot is part of). Long form means that each row in the dataframe only has one value, whereas our data is currently in wide form - there are two values for each individual on each row (V1 & V2). 
We do this using:
```sh
k2_long <- k2 %>% bind_cols(samplelist, k = "k2") %>% 
  pivot_longer(cols = 1:2)
head(k2_long)
```
Lets look at this
- first we call the data (`k2`)
- `%>%` is a way of piping from one command to the next in the package `tidyverse` - similar to the `|` in UNIX
- the function `bind_cols()` is similar to the base R function we used earlier `cbind()` and adds columns together. At the same time we can add a column to specify the value of K
- the function `pivot_longer()` converts it to long form - the number of columns relates to the number of columns in the .Q matrix

![Long_df](../IM/Long_df.png)

> `Exercise six`
>
> Can you modify this code now for your other values of k?

Now we should have a long form data frame object for every cluster with the metadata (samples and regions) attached. 

Using this we will plot the results as a  stacked barplots using `geom_col()`.

We plot in ggplot using the same process as before. First we can make the basic plot:
```sh
admix_plot_k2 <- ggplot(data = k2_long) +
    geom_col(aes(x=sample, y=value, fill=name)) +
    scale_y_continuous(expand = c(0,0))
```
- start with specifying the data
- then we specify the geom and the aesthestics for the bar plot
- the final option (`scale_y_continuous(expand = c(0,0))`) is just to make the bars extend to the bottom of the x-axis

Take a look at the plot, can you tell anything about whether the babirusa from the same region show the same patterns of ancestry? Probably not quite yet as the bars are just in the order R has read the sample column in the metadata. 
We need to split up the columns by the regions, the quick way to do this is with facets. 

In ggplot2, facets are a powerful way to create multiple plots within a single plot. They allow you to group your data by one or more variables and create separate plots for each group. This can be a great way to visualize how your data varies across different groups.

To the base plot add this line: 
```sh
admix_plot_k2 <- admix_plot_k2 + facet_wrap(~region, scales = "free", nrow = 1)
```
- this is the function `facet_wrap()` where we specify we want it to be split by region
- `scales = "free"` stops R trying to plot every sample in every region
- and `nrow = 1` makes sure they are in a line, and not in a square
  
> `Question` - try running it without these options in `facet_wrap()`, what do you see?

Finally to make the names readable and label the y axis correctly we can add: 
```sh
admix_plot_k2 <- admix_plot_k2 + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
 ylab("Admixture proportion")
```

<img src="../IM/ADMIX_K2.png" width="70%" height="70%">

> `Exercise seven`
>
> Save this plot and generate the remaining plots for your clusters?
>
> `Hint` - you can add multiple layer to the object at once
> 
> e.g. `admix_plot_k3 <- ggplot(data = k3_long) +
>                        geom_col(aes(x=sample, y=value, fill=name)) +
>                        scale_y_continuous(expand = c(0,0)) +
>                        facet_wrap(~region, scales = "free", nrow = 1)`
> 

When you have made the plots of all the values of k, we can visualise them ontop of each other using a function called `ggarrange()` in the `ggpubr` package.

You specify the plots you would like to view together, and then the option `ncol = 1` indicates you want them all in one column i.e. ontop of each other.
```sh
admix_all_plots <- ggarrange(admix_plot_k2, admix_plot_k3, admix_plot_k4, admix_plot_k5, ncol = 1)

admix_all_plots
```
If there were any issues install the `ggpubr` package, you can use the `grid.arrange()` a function from the package `gridExtra` which has the same structure
```sh
grid.arrange(admix_plot_k2, admix_plot_k3, admix_plot_k4, admix_plot_k5, ncol = 1)
```
Do not forget to save your script as you go.

<img src="../IM/ADMIX_ALL.png" width="60%" height="60%">

Now we will check to see which value of K the cross validation suggests is best. The easiest way to do this is to return to your terminal and look at the contents of the `cv_errors.txt` file (or you can open this in a text editor). The lower the value the more fitting the k-value is to the data. However, be warned this is only suggestive and should not be overinterpretted.

#### Questions: 
- what is the most fitting value of k? 
- what do you think this analysis tells us about the ancestry of the babirusa?
- are there any individuals that look like they are admixed, i.e. show evidence of multiple ancestries?

### Now that you have the results of all three exercises: 
1. do the analyses agree?
2. what do you think we can conclude overall about the population structure of these babirusa?
3. which region is the most genetically distinct and why might this be?

### 5. Extra exercises
>
> - Can you change the plots to be a better representation of the clusters? `Hint` - look at column we made called `region2` in the metadata. If you use this variable in the plot instead what happens?
>   
> - You can try and the different PC axes along with there specific % contribution of variation for these axes
>   
> - Can you choose your own colour palette for both the PCA and ADMIXTURE? `Hint` - use the `colour=` option in `geom_point()` or the addition layers `scale_colour_manual()`
>  http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
>
> - In the admixture plot can you see that the colours of the clusters are not consistent across the values of k - i.e. V1 in the k2 plot is not the same as V1 in the k3 plot for example
> - There are programs that can be used to fix it - in your browser https://tau.evolseq.net/clumpak/index.html - you will also need to edit the population file to remove the `RD` sample names. > - Ask an instructor if you get to this point
> - Or take a look at the R package http://www.royfrancis.com/pophelper/articles/index.html for some ideas on how to fix this


