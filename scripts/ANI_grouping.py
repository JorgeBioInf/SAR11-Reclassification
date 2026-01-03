"""
ANI_grouping.py
----------------------
Uses ANI information to prepare ConSpeciFix source groups and test groups. 
Given a fastANI output table, splits the genomes into different folders:
    1) n folders with genomes sharing > {ANI_threshold} ANI and with group size > 15 (souce folders)
    2) folder with the remaining genomes (test folder) 

Author: Jorge Marcos Fernández
Date: 2025-11-18
Version: 1.0

Usage:
    python ANI_grouping.py fastANI_results ANI_th

Output:
    - source_genomes_{i}/ folders. i represents as many groups of genomes sharing ANI > ANI_threshold with and size > 15 genomes. 
    - test_genomes/ folder with the rest of the genomes

Dependencies:
    - networkx
    - pandas
    - sys
    - pathlib
    - zipfile
    - os
    - shutil
    - sys
    
Notes:
    - Full directory name found in the names of the genomes in the fastANI output file has to be given in the code 
    for a correct data processing
"""

# To be changed by user
full_dir = "/home/estudiante2/JMF/other_thresholds/SAR11_genomes/"


## LIBRARIES
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
import shutil
import sys
import os
import zipfile

## FUNCTIONS
def copy_file(name, in_folder, out_folder):
    "Copies a file from SAR11 genomes input directory to a target output directory"
    filename = f'{name}.fa'
    in_path = os.path.join(in_folder, filename)
    out_path = os.path.join(out_folder, filename)
    shutil.copyfile(in_path, out_path)
    print(f'Copied {name} --> {out_folder}')

## CHECK ARGUMENTS
args = sys.argv
if len(args) != 3:
    print('Use: ANI_grouping.py fastANI_results_file ANI_threshold')
    sys.exit(1)

df_name = args[1]
colnames = ["Query", "Reference", "ANI", "Bidirectional mappings", "Query fragments"]
try:
    df = pd.read_csv(df_name, sep = '\t', names = colnames, header = None)
except Exception as e:
    print(f'Error reading {df_name}:{e}')
    sys.exit(1)

df["Query"] = df["Query"].str.replace(full_dir, "", regex=False)
df["Query"] = df["Query"].str.replace('.fa', "", regex=False)
df["Reference"] = df["Reference"].str.replace(full_dir, "", regex=False)
df["Reference"] = df["Reference"].str.replace('.fa', "", regex=False)

df.head(10)

# Extract all relations with ANI > threshold
ANI_th = args[2]

try:
    to_group = df[(df["ANI"] >= ANI_th) & (df["Query"] != df["Reference"])]
except:
    print(f"Error: invalid ANI threshold: {ANI_th}")
    sys.exit(1)


uniques = set(to_group["Query"].tolist())
print('Unique sequences found:', len(uniques))
all_uniques = set(df["Query"].tolist())
print('All unique sequences:', len(all_uniques))


### GROUPING

# Create new graph
G=nx.Graph()

# Nodes
for unique in uniques:
    G.add_node(unique)

# Edges
for idx, row in to_group.iterrows():
    node1 = row["Query"]
    node2 = row["Reference"]
    G.add_edge(node1, node2)  


# Connected componets of the graph
components = list(nx.connected_components(G))
components = sorted(components, key = len)

print("Groups found:")
for i, comp in enumerate(components):
    print(f"Group {i+1} ({len(comp)}): {comp}")



## If genome in component with size > 15 --> is a source genome, store in folder source_genomes_X
## Else --> is a test genome, store in folder test_genomes

in_zip = "SAR11_genomes.zip"
in_folder = "SAR11_genomes"

# Uncompress zip with SAR11 genomes 
if not os.path.exists(in_folder):
    with zipfile.ZipFile(in_zip, 'r') as zip_ref:
        zip_ref.extractall(in_folder)
        print(f"{in_zip} extracted to {in_folder}")
        
# Generate test output folder 
test_out_folder = 'test_genomes'
os.makedirs(test_out_folder, exist_ok = True)

# Process source genomes
s = False
source_num = 1
components_list = []

for component in components:
    if len(component) >= 15:  # Source group found
        s = True
        source_out_folder = f'source_genomes_{source_num}'  # Create specific folder
        os.makedirs(source_out_folder, exist_ok = True)
        components_list.append(list(component))
        components_list += 1
        
if not s:
    print('No group of size >= 15 found according to ANI!')
    print('Analysis finished. No output generated')
    exit(1)

## Store genomes
source_genomes = []
for genome in os.listdir(in_folder):
    for i, sublist in enumerate(components_list): # Seach specific source component in which the genome is found
        genome_name = genome.replace('.fa', '')
        if genome_name in sublist:
            source_genomes.append[genome_name]
            copy_file(genome, in_folder, f'source_genomes_{i + 1}') # Copy to specific folder

    if genome not in source_genomes: # If test genome, copy to test folder
        copy_file(genome, in_folder, test_out_folder)

print('\nAnalysis completed! Results are available in folders source_genomes and test_genomes')

