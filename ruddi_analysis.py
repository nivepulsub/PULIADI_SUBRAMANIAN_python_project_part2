#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
import pandas as pd

#using SeqIO module, parse function, parse the fasta file
records = SeqIO.parse("GCF_000287275.1_ASM28727v1_genomic.fna", "fasta")

#creating empty lists to store the required info from iteration
length_of_genome = []
GC_content = []
ATG_forward_count = []
ATG_reverse_count = []

#iterating to get the required info from the fasta file
for record in records: 
    genome = record.seq
    length_of_genome.append(len(genome))
    gc_content_percentage = (gc_fraction(genome))*100
    GC_content.append(gc_content_percentage)
    ATG_forward_count.append(genome.count("ATG"))
    ATG_reverse_count.append(Seq(genome).reverse_complement().count("ATG"))
    
    print(f"{length_of_genome}, {GC_content}, {ATG_forward_count}, {ATG_reverse_count}")
    
#creating dictionary using length_of_genome, length_of_genome, ATG_forward_count, and ATG_reverse_count as keys and their respective values from the empty lists that we created before
dict = {'length_of_genome': length_of_genome, 'GC_content': GC_content, 'ATG_forward_count': ATG_forward_count, 'ATG_reverse_count': ATG_reverse_count}

#creating the dataframe using pandas
df = pd.DataFrame(dict)

#using melt function, converting from wide to long format using pandas
melted_df = pd.melt(df, var_name="Property", value_name="Value")

#writing the dataframe to csv file
melted_df.to_csv('ruddi.csv', index=False)


# In[ ]:




