#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import Entrez
from Bio import SeqIO
import pandas as pd

#Identifying myself to NCBI
Entrez.email = 'nive9756@gmail.com'

#Using Entrez module, efetch function, fetching all the required fasta files
stream = Entrez.efetch(db="protein", id="AGI40145.1,AGJ87295.1,WVV45440.1,WVS05366.1", rettype="fasta", retmode="text")

#parsing all the fasta files using SeqIO module, parse function
records = SeqIO.parse(stream, "fasta")

#creating empty lists to append all the info gather through iteration
accession_list = []
first_10_AA_list = []
Length_list = []
Number_Cs = []

#iterating over to get the id, first 10 aa, length of the protein sequence, and the number of Cysteine for all the fasta files
for record in records:
    accession_list.append(record.id)
    protein_sequence = record.seq
    first_10_aa = protein_sequence[0:10]
    first_10_AA_list.append(first_10_aa)
    protein_sequence_length = len(protein_sequence)
    Length_list.append(protein_sequence_length)
    number_Cs = protein_sequence.count('C')
    Number_Cs.append(number_Cs)
    
    print(f"{record.id}, {first_10_aa}, {protein_sequence_length}, {number_Cs}")
    
#creating dictionary with ID, First_10_AA, Length, and Number_Cs as keys and their respective values from the empty lists that we created before
dict = {'ID': accession_list, 'First_10_AA': first_10_AA_list, 'Length': Length_list, 'Number_Cs': Number_Cs}

#creating dataframe using pandas
df = pd.DataFrame(dict)

#writing dataframe to csv file
df.to_csv('protein_info.csv')


# In[ ]:




