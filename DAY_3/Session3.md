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
Open a new terminal and navigate to the directory where you would like to run the analysis on your local computer

Then run `sftp` - this stands for secure file transfer protocol. We need to tell it where the ssh key is, and which remote server to connect using your username and the IP address. Therefore enter: 
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
- .treefile for the tree
  
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
sftp> bye
```

### 2. Visualising and comparing the trees

### 3. Visualising the PCA analysis
What do you think this tells us about the number of populations

### 4. Visualising the ADMIXTURE analysis
Based on these 










