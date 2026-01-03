"""
GTDB_processer.py
----------------------
From a list of genome files, returns the GTDB taxonomic classification for each identifier. 
Genomes are stored in a separate folder unclassified_gtdb/ if this information is not available for
subsequent analysis with GTDBtk

Author: Jorge Marcos Fernández
Date: 2025-11-30
Version: 1.0

Usage:
    python GTDB_processer.py genomes_directory_path genomes_table_path

Output:
    - GTDB_classification.json dictionary with the GTDB classification of each genome (if available)
    - unclassified_gtdb/ folder with genomes without available GTDB classification

Dependencies:
    - json
    - pandas
    - sys
    - requests
    - os
    - shutil

Notes:
    - Requires directory with genomes FASTAS and table with genomes identifiers (S3 from Free et al)
"""


# LIBRARIES
import json
import requests 
import pandas as pd 
import sys
import os
import shutil


# Check arguments
if len(sys.argv) != 3:
    print('Use: GTDB_processer.ipynb genomes_directory_path genomes_table_path')
    sys.exit(1)

genomes_path = sys.argv[1]
df_path = sys.argv[2]

if not os.path.exists(genomes_path):
    print(f'Error: directory {genomes_path} not found!')
    sys.exit(1)

try:
    df = pd.read_csv(df_path, sep='\t')
except Exception as e:
    print(f'Error reading genomes table: {e}')
    sys.exit(1)

print('Arguments checked successfully')


# Iter the genomes folder 
# If the genome is anotated in GTDB database --> API access 
# If the genome is not anotated --> Store for posterior GTDBtk analysis
print("\nRetrieving GTDB information for genomes")

gtdbtk_genomes = []
gtdb_classification = {}

for file in os.listdir(genomes_path):

    isolate_id = file.split('_')[0]
    name = file.replace('.fa', '')
    refseq_id = df.loc[df["SAG or Isolate ID"] == isolate_id, "RefSeq Assembly (*IMG Genome ID)"].iloc[0]

    url = f'https://gtdb-api.ecogenomic.org/genome/{refseq_id}/taxon-history'
    response = requests.get(url)

    if not response.ok:
        print(f'No information found for {isolate_id}')
        gtdbtk_genomes.append(file)
        continue

    data = response.json()
    
    if len(data) == 0:
        print(f'No information found for {isolate_id}')
        gtdbtk_genomes.append(file)
        continue
    
    print(f'GTDB information obtained for {isolate_id}')
    gtdb_classification[name] = data[0]


## Process data 
# Create folder with unclassified genomes --> unclassified_gtdb
# Create GTDB_classification.json for classified genomes
print('\nGenerating output ...')
cwd = os.getcwd()

out_filename = "GTDB_classification.json"
out_filepath = os.path.join(cwd, out_filename)

with open(out_filepath, 'w') as out_f:
    json.dump(gtdb_classification, out_f, indent=4)

out_dirname = "unclassified_gtdb"
out_dirpath = os.path.join(cwd, out_dirname)
os.makedirs(out_dirpath, exist_ok=True)

for file in gtdbtk_genomes:
    
    full_in_dirpath = os.path.join(genomes_path, file)
    full_out_dirpath = os.path.join(out_dirpath, file)
    
    shutil.copyfile(full_in_dirpath, full_out_dirpath)

print('\nAnalysis completed!')
print(f'You can find GTDB retrieved information in file {out_filename}')
print(f'Unclassified genomes have been copied to folder {out_dirname}')

