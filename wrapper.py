import sys
import os
import argparse

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="PIPELINE PROJECT WRAPPER SCRIPT")

    parser.add_argument("-NCBI", "--NCBI_email",
                        help="Email for accessing NCBI",
                        required=True)

    parser.add_argument("-cwd", "--current_working_directory",
                        help="Current working directory",
                        default=os.getcwd())
    
    parser.add_argument("-R", "--R_directory",
                        help="input path to R file for sleuth",
                        required=True)

    return parser.parse_args(args)

arguments = check_arg(sys.argv[1:])

cwd = arguments.current_working_directory
NCBI_email = arguments.NCBI_email
R_path = arguments.R_directory

########################################################################
from Bio import Entrez
from Bio import SeqIO
import os

index_name = 'index_HCMV.idx' #name of index file
genome_wanted = 'NC_006273.2' #genome we want to access and build transcriptome from 
items_in_directory = os.listdir(cwd) #just a command line command that will list items in your directory for the script to grab 

'''retrieve HCMV genome '''
Entrez.email = NCBI_email #put in your email to use for Entrez 
handle = Entrez.efetch(db="nucleotide", id= genome_wanted, rettype="gb", retmode="text") #then use Efetch to search nucleotides, reference transcriptome of interest, return type, and return mode
record = SeqIO.read(handle, 'genbank') #then read the handle and use genbank format 
handle.close()

########################################################################
'''2 Kallisto'''

from Bio import Entrez
from Bio import SeqIO

'''retrieve HCMV genome '''
Entrez.email = NCBI_email #put in your email to use for Entrez 
handle = Entrez.efetch(db="nucleotide", id= genome_wanted, rettype="gb", retmode="text") #then use Efetch to search nucleotides, reference transcriptome of interest, return type, and return mode
record = SeqIO.read(handle, 'genbank')
handle.close()

'''extract CDS'''
CDS = [] #make a new list to keep our CDS extracted

for items in record.features: #look in record.features
    if items.type == 'CDS': #identify "CDS" type 
        if 'protein_id' in items.qualifiers: #if they have the protein_id in the CDS
            refseq_protein_id = items.qualifiers['protein_id'][0] #then extract the 0th item in protein_id, which is the name
            seq = items.extract(record.seq) #and then extract the sequence SPECIFICALLY for cds
            fasta_format = f">{refseq_protein_id}\n{seq}\n" #write it in the fasta format with > and name and then sequence
            CDS.append(fasta_format) #add each to the CDS list 
            #print(CDS)

with open('reference_transcriptome.fasta', 'w') as f: #then write the fasta file that will contain our fasta sequences 
    for cds_ in CDS: #for each cds...
        f.write(f"{cds_}") #write it into the file, ensures its on a new line every time

total_CDS_found = len(CDS) #find the number of CDS found
with open('PipelineProject.log', 'w') as f: #for the log file, write how many CDS found
    f.write(f"The HCMV genome ({genome_wanted}) has {total_CDS_found} CDS.\n") #format wanted for log file

reference_transcriptome = 'reference_transcriptome.fasta' #keep track of fasta used to index 

import os #run kallisto command automatically to the command line using os.system 
kallisto_command = 'kallisto index -i ' + index_name + ' ' + reference_transcriptome
os.system(kallisto_command)

##################################################################################
'''3'''
with open('PipelineProject.log', 'a') as f: #for the log file, APPEND ('a') putting it now so it doesnt add itself for every sample later.
    f.write(f"sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\n") #format wanted for log file

samples = [] #make an empty list to store your SRR files 
for items in items_in_directory: #for items in your cwd
    if items.startswith('SRR'): #if the files start with SRR (they should because those are the files names i provided)
        samples.append(items) #adds those items to the list to sift 

forward_reads = [] #keeps all the forward reads with _1
reverse_reads = [] #keeps all the reverse reads with _2

for read in samples: #for each item in that list 
    if read.endswith('_1.fastq.gz'): #if the file ends with a underscore 1, then its a forward file and needs to be seperated from the reverse 
        forward_reads.append(read)
        forward_reads.sort() #sort that they will be paired to their corresponding reverse 
    elif read.endswith('_2.fastq.gz'): #same here, if its underscore 2 then it needs to be a reverse file 
        reverse_reads.append(read)
        reverse_reads.sort()

names_list = [] #list that will store the SRR names 
for forward, reverse in zip(forward_reads, reverse_reads): #takes each position into account, pairing up the reads since we sorted them previously so they should be in order
    sample_name = forward.split('_')[0] #takes the first item of the split items, giving the name
    names_list.append(sample_name) #append the name into the file 
    index_file = f'{cwd}/index_HCMV.idx' #this is the name of the index file we created earlier, youll need to access it 

    for names in names_list: #the names correspond to a Donor and Condition, given the ending numbers, we put them into 2 categories Donor 1 or donor 3, and 2dpi and 6dpi
        if names.endswith(('30','33')):
            Donor = 'Donor 1'
        elif names.endswith(('44', '45')):
            Donor = 'Donor 3'
        if names.endswith(('30', '44')):
            condition = '2dpi'
        elif names.endswith(('33', '45')):
            condition = '6dpi'

    output_sample_dir = f'{cwd}/kallisto_results/{sample_name}' #make the output directory so that each paired reads will go into their own file in an output directory 
    os.makedirs(output_sample_dir, exist_ok=True)
    quant_command = "kallisto quant -i " + index_file + " -o " + output_sample_dir + " -b 10 -t 2 " + forward + " " + reverse #quant command requires an input, output, and forward and reverse files 
    os.system(quant_command)

    TPMs = [] #this will keep the TPMs we get from the abundance tsv file output from kallisto 
    with open(f'{output_sample_dir}/abundance.tsv', 'r') as f: #this will go through each abundance tsv individually 
        next(f) #skip the first line since its a header
        for line in f:
            rows = line.strip().split('\t')#strip each line so we dont get other characters, and split by tab since its a tsv
            tpm = float(rows[4]) #the TPMs are the 4th column in the file, idk why i called it row, but im not going to change it now lol 
            TPMs.append(tpm) #then it appends that 4th column item in the TPM 

    with open('PipelineProject.log', 'a') as f: #for the log file, APPEND ('a' the sample, condition, and TPM stats)
        ''' wanted from the tpms '''
        min_TPM = min(TPMs)
        sum_TPM = sum(TPMs)
        mean_TPM = sum(TPMs)/ len(TPMs)
        max_TPM = max(TPMs)

        order_TPMs = sorted(TPMs)
        len_TPMs = len(order_TPMs)
        middle = len_TPMs / 2

        if len_TPMs % 2 == 0: #if the list is even, the median will be the number in the middle 
            half1 = order_TPMs[len_TPMs // 2 - 1]
            half2 = order_TPMs[len_TPMs // 2]
            median = (half1 + half2) / 2 #how the median is calculated for even number lists
        else: #for odd length lists
            median = order_TPMs[len_TPMs // 2] #going to be the middle number index position 

        f.write(f"{sample_name}\t{condition}\t{min_TPM}\t{median}\t{mean_TPM}\t{max_TPM}\n") #format the math results how it wants to the log file 


###############################################################################################
'''4'''
'''rerun code to make sleuth table, and add path to the file.'''

with open('sample_table_final.txt', 'w') as f: #R needs a table of the tpms to make the calulations so we are making one using the code from the kallisto stuff
    f.write(f"sample\tcondition\tmin_tpm\tmed_tpm\tmean_tpm\tmax_tpm\tpath\n") #format wanted for log file


names_list = []
for forward, reverse in zip(forward_reads, reverse_reads): #takes each position into account, pairing up the reads since we sorted them previously so they should be in order
    sample_name = forward.split('_')[0] #takes the first item of the split items, giving the name
    names_list.append(sample_name)
    index_file = f'{cwd}/index_HCMV.idx'
    path = f'{cwd}/kallisto_results/{sample_name}' #path for R script, it needs a path to access the kallisto output files

    for names in names_list: #repeat but cant hurt 
        if names.endswith(('30','33')):
            Donor = 'Donor 1'
        elif names.endswith(('44', '45')):
            Donor = 'Donor 3'
        if names.endswith(('30', '44')):
            condition = '2dpi'
        elif names.endswith(('33', '45')):
            condition = '6dpi'

    TPMs = []
    with open(f'{output_sample_dir}/abundance.tsv', 'r') as f: #access tpms, see above in #3 for more thorough comments 
        next(f) #skip the first line
        for line in f:
            rows = line.strip().split('\t')#strip each line so we dont get other characters, and split by tab since its a tsv
            tpm = float(rows[4])
            TPMs.append(tpm)

    with open('sample_table_final.txt', 'a') as f: #for the log file, APPEND ('a' the sample, condition, and TPM stats)

        min_TPM = min(TPMs)
        sum_TPM = sum(TPMs)
        mean_TPM = sum(TPMs)/ len(TPMs)
        max_TPM = max(TPMs)

        order_TPMs = sorted(TPMs)
        len_TPMs = len(order_TPMs)
        middle = len_TPMs / 2

        if len_TPMs % 2 == 0: #if the list is even, the median will be the number in the middle 
            half1 = order_TPMs[len_TPMs // 2 - 1]
            half2 = order_TPMs[len_TPMs // 2]
            median = (half1 + half2) / 2 #how the median is calculated for even number lists
        else: #for odd length lists
            median = order_TPMs[len_TPMs // 2] #going to be the middle number index position 

        f.write(f"{sample_name}\t{condition}\t{min_TPM}\t{median}\t{mean_TPM}\t{max_TPM}\t{path}\n")


r_sleuth_command = f'Rscript {R_path}' #use Rscript command to run the rscript in a python script 
os.system(r_sleuth_command)


#################################################################################################
'''5'''

'''retrieve HCMV genome '''
Entrez.email = NCBI_email #put in your email to use for Entrez 
handle = Entrez.efetch(db="nucleotide", id= genome_wanted, rettype="gb", retmode="text") #then use Efetch to search nucleotides, reference transcriptome of interest, return type, and return mode
record = SeqIO.read(handle, 'genbank')
handle.close()

genome_sequence = record.seq #make a genome file, but we just want it in a fasta with the name and the sequence 
with open('HCMV_genome_index.fasta', 'w') as f:
    f.write(f'> {genome_wanted}\n{genome_sequence}')

bowtie_index_command = f'bowtie2-build {cwd}/HCMV_genome_index.fasta HCMV_index_bowtie' #bowtie command to make an index genome 
os.system(bowtie_index_command)

'''bowtie2 assembly'''
samples = []
for items in items_in_directory:
    if items.startswith('SRR'):
        samples.append(items)
        #print(samples)

forward_reads = [] #keeps all the forward reads with _1
reverse_reads = [] #keeps all the reverse reads with _2

for read in samples:
    if read.endswith('_1.fastq.gz'):
        forward_reads.append(read)
        forward_reads.sort()
    elif read.endswith('_2.fastq.gz'):
        reverse_reads.append(read)
        reverse_reads.sort()


names_list = []
for forward, reverse in zip(forward_reads, reverse_reads): #takes each position into account, pairing up the reads since we sorted them previously so they should be in order
    sample_name = forward.split('_')[0] #takes the first item of the split items, giving the name
    names_list.append(sample_name)
    
    for names in names_list:
        if names.endswith(('30','33')):
            Donor = 'Donor 1'
        elif names.endswith(('44', '45')):
            Donor = 'Donor 3'
        if names.endswith(('30', '44')):
            condition = '2dpi'
        elif names.endswith(('33', '45')):
            condition = '6dpi'

    output_sample_dir = f'{cwd}/kallisto_results/{sample_name}' #make the output directory so that each paired reads will go into their own file in an output directory 

    forward_unzipped = f'{sample_name}_1.fastq' #this is what we will call the unzipped fastq because they need to be unzipped for bowtie
    reverse_unzipped = f'{sample_name}_2.fastq'

    os.system(f'gunzip -c {forward} > {forward_unzipped}') #then use gunzup to unzip gzip files. -c will keep them from being overwritten 
    os.system(f'gunzip -c {reverse} > {reverse_unzipped}')

    os.system(f'wc -l {forward_unzipped} > reads_before.txt') #write wc -l results to a file for before filtering results
    with open('reads_before.txt', 'r') as before_file: #read it cause it should be read number and then the sample 
        count_before = int(before_file.read().split()[0]) // 4 #floor divide by 4 because fastq is every 4 lines (4 lines long)

    bowtie2_command_assembly = f"nohup bowtie2 --quiet -x HCMV_index_bowtie -1 {forward_unzipped} -2 {reverse_unzipped} --al-conc {cwd}/mapped_{sample_name}.fastq -S HCMVmap.sam "
    os.system(bowtie2_command_assembly)

    os.system(f'wc -l mapped_{sample_name}.1.fastq > forward_reads_after.txt') #write wc -l results to a file for after filtering results, only need to check one file (either forward or reverse)
    with open('forward_reads_after.txt', 'r') as after_file: #read it cause it should be read number and then the sample 
        count_after = int(after_file.read().split()[0]) // 4 #floor divide by 4 because fastq is every 4 lines (4 lines long)

    with open('PipelineProject.log', 'a') as f: #for the log file
        f.write(f"{Donor} ({condition}) had {count_before} read pairs before Bowtie2 filtering and {count_after} read pairs after. \n") #format the sentence for each sample into the log file 

#################################################################################################
'''6 SPAdes'''
items_in_directory = os.listdir(cwd)
SPA_samples = [] #spades samples will be kept here
for items in items_in_directory:
    if items.startswith('mapped_'): #since we renamed them in the bowtie section as mapped_, look for those in the cwd 
        SPA_samples.append(items)

Donor_1_SPA_samples = [] #for spades we are separating Donor 1 assembly from donor 3 assembly so we make 2 lists 
Donor_3_SPA_samples = []

for samples in SPA_samples:
    if samples.endswith(('30.1.fastq','33.1.fastq','30.2.fastq','33.2.fastq')): #sample 30 and 33 belong to donor 1
        Donor = 'Donor 1'
        Donor_1_SPA_samples.append(samples)
    elif samples.endswith(('44.1.fastq', '45.1.fastq', '44.2.fastq', '45.2.fastq')): #sample 44 and 45 belong to donor 3, so we add them to the appropriate lists 
        Donor = 'Donor 3'
        Donor_3_SPA_samples.append(samples)


D1_SPA_forward_reads = [] #keeps all the mapped forward reads with .1 
D1_SPA_reverse_reads = [] #keeps all the mapped reverse reads with .2

for read in Donor_1_SPA_samples: #now we seperate the sampels in the donor 1 by forward and reverse to use in the spades command 
    if read.endswith('.1.fastq'):
        D1_SPA_forward_reads.append(read)
        D1_SPA_forward_reads.sort()
    elif read.endswith('.2.fastq'):
        D1_SPA_reverse_reads.append(read)
        D1_SPA_reverse_reads.sort()

'''repeat same code but for the other donor'''

D3_SPA_forward_reads = [] #keeps all the mapped forward reads with .1
D3_SPA_reverse_reads = [] #keeps all the mapped reverse reads with .2

for read in Donor_3_SPA_samples: #now we seperate the sampels in the donor 3 by forward and reverse to use in the spades command 
    if read.endswith('.1.fastq'):
        D3_SPA_forward_reads.append(read)
        D3_SPA_forward_reads.sort()
    elif read.endswith('.2.fastq'):
        D3_SPA_reverse_reads.append(read)
        D3_SPA_reverse_reads.sort()

with open('PipelineProject.log', 'a') as f: #now we run that and add the command used to the log file 
    for forw1, rev1 in zip(D1_SPA_forward_reads, D1_SPA_reverse_reads): #for the forward and reverse reads for donor 1 in their lists as a pair...
        spades_command_D1 = f'spades.py -k 77 -t 2 --only-assembler -1 {forw1} -2 {rev1} -o {cwd}/SPAdes_assembly/Donor1'
        os.system(spades_command_D1)
        f.write(f"{spades_command_D1}\n") #should be 2 commands total 


    for forw3, rev3 in zip(D3_SPA_forward_reads, D3_SPA_reverse_reads):
        spades_command_D3 = f'spades.py -k 77 -t 2 --only-assembler -1 {forw3} -2 {rev3} -o {cwd}/SPAdes_assembly/Donor3' #for the forward and reverse reads for donor 3 in their lists as a pair...
        os.system(spades_command_D3)
        f.write(f"{spades_command_D3}\n") #should be 2 commands total 

        #4 spades commands total should be added to log file 

###########################################
'''7 BLAST'''

longest_contig = 0 #the longest contig length will be incremented so we start with 0 to initalize 
with open(f'{cwd}/SPAdes_assembly/Donor1/contigs.fasta', 'r') as file: #then it looks though all the contigs made by spades for donor 1
    for record in SeqIO.parse(file, 'fasta'): #since they are a fasta, parse fasta 
        if len(record.seq) > longest_contig: #and if it finds a sequence longer than the current longest contig
            longest_contig_D1 = record.seq #then that sequence becomes the longest contig for D1
            longest_contig = len(record.seq) #and its length gets stored and we restart 

            with open('Donor1_contig.fasta', 'w') as f: #then we make the longest contig a fasta file so we can use it in blast
                f.write(f"> Donor1\n{longest_contig_D1}\n")

## the same as seen above is done but for donor 2. we will have 2 fasta files in total with 2 longest contigs one for each donor 
longest_contig = 0 
with open(f'{cwd}/SPAdes_assembly/Donor3/contigs.fasta', 'r') as file:
    for record in SeqIO.parse(file, 'fasta'):
        if len(record.seq) > longest_contig:
            longest_contig_D3 = record.seq
            longest_contig = len(record.seq)

            with open('Donor3_contig.fasta', 'w') as f:
                f.write(f"> Donor3\n{longest_contig_D3}\n")



#making Betaherpesvirinae nucleotide database 
Betaherpesvirinae_DB_command = 'datasets download virus genome taxon Betaherpesvirinae --include genome' #making a database, but first we have to download all the Betaherpesvirinae data!!
os.system(Betaherpesvirinae_DB_command)

unzip_ncbi = f'unzip {cwd}/ncbi_dataset.zip' #all the info for Betaherpesvirinae is stored in the ncbi zip so we have to unzip it to use it 
os.system(unzip_ncbi)

make_DB = f'makeblastdb -in {cwd}/ncbi_dataset/data/genomic.fna -out Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl' #then we make the database using the data we just downloaded for Betaherpesvirinae
os.system(make_DB)

#then we use the blast command, blastn for nucleotide each longest contig for each donor. we only want the top 10 so i used -max_target_seqs to limit to 10, max_hsps since we want the BEST alignments (HSP) and then the output file and the output format which was asked for in the prompt 
BLAST_D1 = f'blastn -query {cwd}/Donor1_contig.fasta -db Betaherpesvirinae -max_target_seqs 10 -max_hsps 1 -out Donor1_results.csv -outfmt "10 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"'
BLAST_D3 = f'blastn -query {cwd}/Donor3_contig.fasta -db Betaherpesvirinae -max_target_seqs 10 -max_hsps 1 -out Donor3_results.csv -outfmt "10 qseqid sacc pident length qstart qend sstart send bitscore evalue stitle"'

os.system(BLAST_D1)
os.system(BLAST_D3)

with open('PipelineProject.log', 'a') as f: #then we are going to append the blast results to the log file
    with open(f'{cwd}/Donor1_results.csv', 'r') as csv: #the blast results are in the CSV
        results_D1 = csv.readlines() #read all the lines 
        header_inlog = False #and Donor 1/ donor 3 so going to be the header, but i have to put false so it doesnt rewrite it every time

        for result in results_D1:  #for the items in the CSV
            each= result.strip().split('\n') #each result  has a new line character so we split them using that

            for item in each: #and then for each result
                item = item.split(',') #split up the necessary components using split but by the comma 
                name = item[0]
                sacc = item[1]
                pident = item[2]
                length = item[3]
                qstart = item[4]
                qend = item[5]
                sstart = item[6]
                send = item[7]
                bitscore = item[8]
                evalue = item[9]
                stitle = item[10]

                if not header_inlog: #and if the header isnt there yet, add it 
                    f.write(f"{name}\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                    header_inlog = True #and then mark it true so it doesnt add it again 
                
                f.write(f"{sacc}\t{pident}\t{length}\t{qstart}\t{qend}\t{sstart}\t{send}\t{bitscore}\t{evalue}\t{stitle}\n") #and then write all the results in a tab delimited order

#do the same as above but just with donor 3, see above for more detailed comments 
    with open(f'{cwd}/Donor3_results.csv', 'r') as csv:
        results_D3 = csv.readlines()
        header_inlog = False

        for result in results_D3:
            each= result.strip().split('\n')

            for item in each:
                item = item.split(',')
                name = item[0]
                sacc = item[1]
                pident = item[2]
                length = item[3]
                qstart = item[4]
                qend = item[5]
                sstart = item[6]
                send = item[7]
                bitscore = item[8]
                evalue = item[9]
                stitle = item[10]

                if not header_inlog:
                    f.write(f"{name}\nsacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                    header_inlog = True
                
                f.write(f"{sacc}\t{pident}\t{length}\t{qstart}\t{qend}\t{sstart}\t{send}\t{bitscore}\t{evalue}\t{stitle}\n")


#the end! 