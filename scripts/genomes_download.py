"""
genomes_download.py
----------------------
Retrieves FASTA sequence for all SAR11 genomes registered in the Supplementary Table S3 from by Free et. al (2024),
considering specific completeness and contamination thresholds. 

Author: Jorge Marcos Fernández
Date: 2025-11-07
Version: 1.0

Usage:
    python genomes_download.py

Output:
    - SAR11_genomes.zip with FASTA sequences of filtered genomes
    - IMG_Genome_IDs.txt with identifiers of the genomes if their FASTA sequences were not found in RefSeq

Dependencies:
    - os
    - pandas
    - requests
    - shutil
    - zipfile

Notes:
    - Accepts no arguments
    - Contamination and completeness thresholds must be changed directly in the code
    - Requires SAR11_genomes_1.zip with all the genomes from the reference article, which can be directly dowloaded from 
      https://figshare.com/articles/dataset/New_SAR11_isolate_genomes_from_the_tropical_Pacific_Ocean/28087454/1
    - Requires Genomes_table.txt with the supplementary table S3 from the reference article, which can be downloaded from
      https://figshare.com/articles/dataset/Supplementary_tables_for_SAR11_genomes_from_the_tropical_Pacific_study_Freel_et_al_2024_/28087490/1?file=51364793
"""

# Completeness and contamination thresholds - might be changed by the user
comp_th = 90
cont_th = 5


## LIBRARIES
import pandas as pd 
import requests
import zipfile
import os
import shutil


## FUNCTIONS
def data_filtering(data, comp_th, cont_th):
    """Filters a dataset accoriding to completeness and contamination thresholds"""

    data = data.copy()
    
    group_list = list(data["Subgroup"])
    groups = list(set(data["Subgroup"].dropna()))

    # Transform into numeric
    data["Completeness"] = data["Completeness"].astype(str).str.replace(',', '.')
    data["Contamination"] = data["Contamination"].astype(str).str.replace(',', '.')

    data["Completeness"] = pd.to_numeric(data["Completeness"], errors='coerce')
    data["Contamination"] = pd.to_numeric(data["Contamination"], errors='coerce')

    # Apply filters
    sset = data[(data["Completeness"] >= comp_th) & (data["Contamination"] <= cont_th)]
    group_sset_list = list(sset["Subgroup"])

    print('Total maintained:', len(group_sset_list),'\n')
    
    # Calculate metrics
    for g in groups:
        count = group_list.count(g)
        sset_count = group_sset_list.count(g)

        fraction = f'{sset_count}/{count}'
        percentage = round(sset_count / count * 100,2)

        print(f'{fraction} ({percentage}%) samples are maintained for group {g}!')
    
    return sset


## MAIN PROGRAM 

### Read RefSeq data
df = pd.read_csv("Genomes_table.txt", sep = '\t')
df = df[(~(df["Category"] == "Outgroup")) & (~(df["Category"] == "This study"))]
print(df.shape)
df.head()


# Filter data
print('\nREFSEQ DATA\n')
filt_data = data_filtering(df, comp_th, cont_th)


# Get RefSeqs 
ids = filt_data['RefSeq Assembly (*IMG Genome ID)']
ids = ids[ids.notna()].tolist()
print(f'{len(ids)} genomes found!\n')


# Filter RefSeq IDs and IMG Genome ID
IMG_IDs = []  

for ID in ids:
    if ID.endswith('*'):
        IMG_IDs.append(ID)
        ids.remove(ID)

filename = 'IMG_Genome_IDs.txt'

if IMG_IDs:
    with open(filename, 'w') as f:
        for i in IMG_IDs:
            f.write(f'{i}\n')

    print(f'{len(IMG_IDs)} IMG Genome IDs found! They must be manually downloaded')
    print('They can be found in IMG_Genome_IDs.txt\n')


# Output folder
path = os.path.join(os.getcwd(), 'SAR11_genomes/')

if not os.path.exists(path):
    os.makedirs(path, exist_ok=True)


# Retrieve sequences
RefSeq = ids
for acc in RefSeq:

    filename = f'{acc}.zip'
    original_path = os.path.join(os.getcwd(), acc)
    fasta_subpath = f'{acc}/ncbi_dataset/data/{acc}/'
    fasta_path = os.path.join(os.getcwd(), fasta_subpath)
    
    # Access genome files
    url = (
        f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download'
        f'?include_annotation_type=GENOME_FASTA'
    )

    response = requests.get(url)

    if not response.ok:
        print(f'There was an error when searching {acc} in RefSeq')
    
    # Write zip 
    with open(filename, 'wb') as file:
        file.write(response.content)
          
    # Unzip information
    try:
        with zipfile.ZipFile(filename, "r") as zip_ref:
            zip_ref.extractall(acc)
            print(f"Extracted: {acc}/")        
        os.remove(filename)
    
    except:
        print(f'Error reading {filename}')
        print(f'Skipping ...')
        continue

    if not os.path.exists(fasta_path):
        print(f'Error: no fasta file found for {acc}')
        continue
        
    # Move fasta to common directory
    for f in os.listdir(fasta_path):
        if f.endswith('.fa') or f.endswith('.fasta') or f.endswith('.fna'):
            SAG = filt_data.loc[filt_data["RefSeq Assembly (*IMG Genome ID)"] == acc, "SAG or Isolate ID"].iloc[0]
            group = filt_data.loc[filt_data["RefSeq Assembly (*IMG Genome ID)"] == acc, "Subgroup"].iloc[0]
            new_filename = f'{SAG}_{group}.fa'
            new_filepath = os.path.join(path, new_filename)
            shutil.move(os.path.join(fasta_path, f), new_filepath)
            print(f'Succesfully moved: {f} --> {path}\n')
    
    # Delete info
    shutil.rmtree(original_path)    


# ### Article data
df2 = pd.read_csv("Genomes_table.txt", sep = '\t')
sset2 = df2[df2["Category"] == "This study"]
study_data_filt = data_filtering(sset2, comp_th, cont_th)


# Delete all study genomes that do not follow criteria
isolate_IDs = sset2["SAG or Isolate ID"].tolist()
filt_isolate_IDs = study_data_filt["SAG or Isolate ID"].tolist()
to_delete_IDs = list(set(isolate_IDs)-set(filt_isolate_IDs))

# Get directory with SAR11_Genomes_1.zip
study_genomes_dirname = 'SAR11_Genomes_1.zip'
study_genomes_path = os.path.join(os.getcwd(), study_genomes_dirname)

# Unzip folder
with zipfile.ZipFile(study_genomes_path, "r") as zip_ref:
        zip_ref.extractall(path)

# Remove specific fasta files
for name in to_delete_IDs:
    
    filename = f'{name}.fa'
    fasta_path = os.path.join(path, filename)
    
    os.remove(fasta_path)

# Rename remaining files to add group info
for name in filt_isolate_IDs:
    
    filename = f'{name}.fa'
    fasta_path = os.path.join(path, filename)
    
    group = study_data_filt.loc[study_data_filt["SAG or Isolate ID"] == name, "Subgroup"].iloc[0]
    new_filename = f'{name}_{group}.fa'
    new_filepath = os.path.join(path, new_filename)
    
    shutil.move(fasta_path, new_filepath)


# In[ ]:


output_name = "SAR11_genomes"
shutil.make_archive(output_name, 'zip', path)
print(f"\nFolder: {output_name}.zip created!")

shutil.rmtree(path) 

