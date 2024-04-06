#!/usr/bin/env python
# coding: utf-8

# In[3]:


from Bio import Entrez
from Bio import SeqIO
import pandas as pd

#Identifying myself to NCBI
Entrez.email = 'nive9756@gmail.com'

#Using Entrez module, efetch function, fetching all the required genbank files
stream = Entrez.efetch(db="nucleotide", id="NZ_CALPCP010000001.1,NZ_CALPCY010000130.1,NZ_BHVZ01000001.1,NZ_SRYA01000017.1,NZ_CAJTFZ010000019.1", rettype="gb", retmode="text")

#parsing using SeqIO module, parse function 
records = SeqIO.parse(stream, "gb")

#creating ampty list to append the info from all the genbank files
accession_list = []
family_list = []
genus_list = []
species_list = []
num_features_list = []
source_list = []

#iterate to get the required features and annotations
for record in records:
    accession_list.append(record.id)
    
    organism = record.annotations['organism']
    [genus, species] = organism.split()
    genus_list.append(genus)
    species_list.append(species)
    
    taxonomy_list = record.annotations['taxonomy']
    family = taxonomy_list[-2] #-1 has genus, -2 is family
    family_list.append(family)
    
    num_features_list.append(len(record.features))
    source_list.append(record.annotations['source'])
    
    print(f"{record.id}, {family}, {genus}, {species}, {len(record.features)}, {record.annotations['source']}")

#creating dictionary with Accession, Family, Genus, Species, Num_Faeatures, and Source as keys and their respective values from the list that we created before 
dict = {'Accession': accession_list, 'Family': family_list, 'Genus': genus_list, 'Species': species_list, 'Num Features': num_features_list, 'Source': source_list}

#creating dataframe using pandas
df = pd.DataFrame(dict)

#writing the dataframe to csv file 
df.to_csv('genbank_parse.csv')


# In[ ]:





# In[ ]:




