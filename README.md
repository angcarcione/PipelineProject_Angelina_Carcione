# Pipeline_Project
NCBI_email = input('Input your NCBI email:')
reference_transcriptome = input('Input accession to build transcriptome:') #'NC_006273.2'

from Bio import Entrez
from Bio import SeqIO

Entrez.email = NCBI_email #put in your email to use for Entrez 
handle = Entrez.efetch(db=nucleotide, id= reference_transcriptome, rettype=fasta, retmode=text) #then use Efetch to search nucleotides, reference transcriptome of interest, return type, and return mode
record = SeqIO.read(handle, 'fasta')
handle.close()
#print(handle)
