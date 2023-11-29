########################
### Installing conda ###
########################
cd /home/
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6 libxml2-dev libcurl4-openssl-dev libfontconfig1-dev libssl-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev

sudo curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
sudo bash /home/Anaconda3-2023.09-0-Linux-x86_64.sh
# say yes to licence agreement, enter /home/anaconda3 as path, say no to initialization


#############################
### Creating conda groups ###
#############################
sudo groupadd anaconda
sudo chgrp -R anaconda /home/anaconda3
sudo chown -R root:anaconda ./anaconda3/
chmod 775 -R /home/anaconda3
sudo adduser ubuntu anaconda
source /home/anaconda3/bin/activate
conda init

########################################
### Adding users to the conda group ####
########################################
while read -r line
sudo adduser $line anaconda
done< /home/ubuntu/Participants/users_list.txt

## note that user need to be created manually first 
## using the adduser command and fill the correct password
## to use conda each needs to run the source command above
## which is part of the tutorial ;)


#########################################
### Adding instructors to the sudoers ###
#########################################
sudo usermod -aG sudo sabhrina1
sudo usermod -aG sudo rosie1
sudo usermod -aG sudo peter1
sudo usermod -aG sudo rasmus1
sudo usermod -aG sudo anubhab1
sudo usermod -aG sudo laurent1

#####################################################
### Creating shared folder and changing ownership ###
#####################################################
cd /home/
sudo mkdir DATA
cd DATA
sudo mkdir Day_1
sudo mkdir Day_2
sudo mkdir Day_3_a
sudo mkdir Day_3_b
sudo mkdir Day_4
sudo mkdir Day_5_s1
sudo mkdir Day_5_s2
sudo mkdir Day_5_s3

sudo chown -R ubuntu:anaconda Day_1/
sudo chown -R ubuntu:anaconda Day_2/
sudo chown -R anubhab1:anaconda Day_3_a/
sudo chown -R rosie1:anaconda Day_3_b/
sudo chown -R sabhrina1:anaconda Day_4/
sudo chown -R ubuntu:anaconda Day_5_s1/
sudo chown -R rasmus1:anaconda Day_5_s2/
sudo chown -R rosie1:anaconda Day_5_s3/


###########################
### Creating conda envs ###
###########################
conda create --name Day_1
conda create --name Day_2
conda create --name Day_3y
conda create --name Day_4
conda create --name Day_5_s1
conda create --name Day_5_s2
conda create --name Day_5_s3

#############
### Day_1 ###
#############
conda activate Day_1
conda install -c bioconda plink
conda install -c bioconda bedtools
conda deactivate

#############
### Day_2 ###
#############
conda activate Day_2
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda adapterremoval
conda install -c bioconda plink
conda install -c bioconda freebayes
conda deactivate

#############
### Day 3 ###
#############
conda activate Day_3
conda install -c conda-forge bcftools=1.18
conda install -c bioconda plink
conda install -c bioconda iqtree
conda install -c bioconda vcftools
conda install -c bioconda eigensoft
conda install -c bioconda admixture
conda install -c bioconda admixtools
conda deactivate

#############
### Day 4 ###
#############
conda activate Day_4
conda install -c conda-forge bcftools=1.18
conda install -c bioconda plink
conda install -c bioconda vcftools
conda install -c genomedk psmc
conda install -c bioconda bedtools
conda install texlive-core
conda install -c conda-forge gnuplot
conda install -c conda-forge ghostscript
conda deactivate

#############
### Day 5 ###
#############
conda activate Day_
conda install -c bioconda admixtools
conda install -c bioconda plink
conda deactivate


#########################################################
### Installing python2 and python3 modules systemwide ###
#########################################################
sudo apt-add-repository universe
sudo apt update
sudo apt install python2-minimal
sudo apt install python3-pip
sudo pip install numpy
sudo pip install pandas

##############################################
### Installing R and R-packages systemwide ###
##############################################
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt update
sudo apt install r-base
sudo -i R
# now we can install all required packages
