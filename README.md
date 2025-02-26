# Pipeline_Project (Angelina Carcione)

In this repository you will find the necessary files and dependencies for running the wrapper script for this pipeline. 
There is also a PipelineProject.log file that shows the results of the script using the *original* data.

The pipeline is as follows: 
- download SRA files
- (KALLISTO) build a transcriptome index for HCMV (NCBI accession NC_006273.2) through CDS extraction
- (KALLISTO) quantify the TPM of each CDS
- (SLEUTH) find differentially expressed genes between conditions (2dpi/6dpi)
- (BOWTIE2) index genome  HCMV (NCBI accession NC_006273.2) and assemble reads
- (SPAdes) using reads that assembled, assemble based on Donor (Donor1/Donor3)
- (BLAST) using the assembly from each donor, BLAST longest contig and identify best alignments

######
# **How to run the wrapper.py**

To run the python script with the tools, use this command line function with REQUIRED flags.

Required Flags: 
- cwd = path to your current working directory (this will be the directory with the repository, you must CD into said directory before running script)
- NCBI_email = email used to access NCBI 
- R_path = path to R script included in repository (Sleuth.R)

**cloning the repo to run the script**
```
git config --global user.name '#YOUR NAME'
git config --global user.email '#YOUREMAIL@gmail.com'

# as a check 
git config --global --list

# should return 
user.name= YOUR NAME
user.email= YOUREMAIL@gmail.com
```
COPY THIS EXACTLY: 
```
git clone https://github.com/angcarcione/PipelineProject_Angelina_Carcione.git
cd PipelineProject_Angelina_Carcione
ls
cat README.md
git remote show origin
```
**YOU MUST CD INTO THE REPO ONCE CLONING**

## COMMAND TO RUN ONCE YOU HAVE CLONED REPO (~ will be your home directory)
**python ~/PipelineProject_Angelina_Carcione/wrapper.py -NCBI YOUREMAIL@luc.edu -cwd ~/PipelineProject_Angelina_Carcione -R ~/PipelineProject_Angelina_Carcione/Sleuth.R**


what mine looks like (DONT USE THIS): python /home/2025/acarcione/PipelineProject_Angelina_Carcione/wrapper.py -NCBI acarcione@luc.edu -cwd /home/2025/acarcione/PipelineProject_Angelina_Carcione -R /home/2025/acarcione/PipelineProject_Angelina_Carcione/Sleuth.R 

*PLEASE NOTE*

Once you run, you will not need to input anything after running the above command EXCEPT for when it comes to BLAST and Downloading: ncbi_dataset.zip.

It will ask you if it is okay to replace the READme.md. You can say whatever you deem fit, I say yes (y) as it will not affect the results. 

######
# How the real SRA files were fetched

If you are interested in running the script with the real SRA files (not the sample data), this is how to retrieve it: 

'''to get the SRR ids, go to https://www.ebi.ac.uk/ena/browser/search or https://sra-explorer.info/ 
and type in the experiment IDs (start with SRX) and download the script to the forward and reverse fastq files'''
```
import os
wget_command = "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/000/SRR5660030/SRR5660030_1.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/000/SRR5660030/SRR5660030_2.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/003/SRR5660033/SRR5660033_1.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/003/SRR5660033/SRR5660033_2.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/004/SRR5660044/SRR5660044_1.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/004/SRR5660044/SRR5660044_2.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/005/SRR5660045/SRR5660045_1.fastq.gz \
     -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR566/005/SRR5660045/SRR5660045_2.fastq.gz "

os.system(wget_command)
```
running the above script in python will download the gzipped files. 

######
# How the sample data was made

Sample data is a truncated version of the SRA data seen above. 

I used command: 
```
head -n 40000 current_data.fastq > sampledata.fastq
```
i then used 'gzip' to zip them for storage purposes as well as how the files will be run in the script.

######
# Dependencies needed and how to download them 

If you are not on a server that has the tools already you will need to download:
- **Kallisto**
https://pachterlab.github.io/kallisto/download (other download options on website) 
with conda/bioconda: conda install kallisto

- **Bowtie2**
https://github.com/BenLangmead/bowtie2 (other download options on website) 
with conda/bioconda: conda install bowtie2

- **Sleuth** (library(sleuth) in R, already included in Sleuth.R script)
it is already included in the R script provided but youll need:
library(sleuth) and library(dplyr)

- **SPAdes**
https://ablab.github.io/spades/installation.html
run as a spades.py script
for linux: https://ablab.github.io/spades/installation.html#downloading-spades-linux-binaries

- **NBCI's BLAST+**
go to the following website and download the appropriate version for your computer: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
