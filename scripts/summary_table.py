"""
summary_table.py
----------------------
Merges information from all classification approaches into a summary table.
For each genome, provides the following information:
    - Accession number
    - Isolate ID
    - Clade
    - % Contamination
    - % Completeness
    - ANI specie
    - GTDB specie
    - ConSpeciFix specie
    - PopCOGenT specie
    - Proposed specie

Author: Jorge Marcos Fernández
Date: 2025-12-20
Version: 1.0

Usage:
    python Summary_table.py SAR11_genomes_list.txt genomes_table supplementary_data_table ANI_table GTDB_classification.json 
    PopCOGenT_table CSF_results.json number_of_CSF_sources

Output:
    - genomes_classification.tsv table

Dependencies:
    - networkx
    - pandas
    - sys
    - json
    
Notes:
    - Requires all classification workflow to be previously executed
"""


# PACKAGES
import pandas as pd  
import networkx as nx 
import sys
import json

# CHECK ARGUMENTS
num_args = len(sys.argv)

if num_args != 9:
    print(f'Error: 8 arguments needed; {num_args - 1} given')
    print(f'Use: python Summary_table.py <SAR11_genomes_list.txt> <genomes_table> <supplementary_data_table> <ANI_table> <GTDB_classification.json> <PopCOGenT_table> <CSF_results.json> <number_of_CSF_sources>')
    sys.exit(1)

# READ DATA

# All genome names 
all_genomes_file = sys.argv[1]

try: 
    with open(all_genomes_file, 'r') as file:
        all_genomes = [line.strip() for line in file]
except Exception as e:
    print(f'Error reading genomes name: {e}')
    sys.exit(1)


# Genomes info
genomes_table = sys.argv[2]

try:
    gen_df = pd.read_csv(genomes_table, sep = '\t')
except Exception as e:
    print(f'Error reading genomes table: {e}')
    sys.exit(1)

# Supplementary info from article
suplementary_table = sys.argv[3]

try:
    sup_df = pd.read_csv(genomes_table, sep = '\t')
except Exception as e:
    print(f'Error reading suplementary table: {e}')
    sys.exit(1)


# ANI information
ANI_table = sys.argv[4]

try:
    colnames = ["Query", "Reference", "ANI", "Bidirectional mappings", "Query fragments"]
    df = pd.read_csv(ANI_table, sep = '\t', names = colnames, header = None)
except Exception as e:
    print(f'Error loading ANI table: {e}')
    sys.exit(1)

# GTDB information
GTDB_json = sys.argv[5]

try:
    with open(GTDB_json, 'r') as file:
        GTDB_data = json.load(file)
except Exception as e:
    print(f'Error reading GTDB information: {e}')
    sys.exit(1)

# PopCOGenT information
POP_data = sys.argv[6]

try:
    pop_df = pd.read_csv(POP_data, sep = '\t')
except Exception as e:
    print(f'Error reading PopCOGenT table: {e}')
    sys.exit(1)

# ConSpeciFix information
CSF_json = sys.argv[7]

try:
    with open(CSF_json, 'r') as file:
        CSF_data = json.load(CSF_json)
except Exception as e:
    print(f'Error reading ConSpeciFix information: {e}')
    sys.exit(1)

sources_num = int(sys.argv[8])  # Number of source genomes

print('All data read successfully!\n')


### ANI grouping
full_dir = "/home/estudiante2/JMF/other_thresholds/SAR11_genomes/"  # Might be changed depending on ANI data

df["Query"] = df["Query"].str.replace(full_dir, "", regex=False)
df["Query"] = df["Query"].str.replace('.fa', "", regex=False)
df["Reference"] = df["Reference"].str.replace(full_dir, "", regex=False)
df["Reference"] = df["Reference"].str.replace('.fa', "", regex=False)

# Extract all relations with > 95% ANI
to_group = df[(df["ANI"] >= 95) & (df["Query"] != df["Reference"])]
uniques = set(to_group["Query"].tolist())

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

components = sorted(nx.connected_components(G), key=len, reverse=True)
sources = components[0:sources_num]
groups = components[sources_num:len(components)]


# ==========================

### MAIN 
accs = []
isolate_ids = []
clades = []
comps = []
conts = []
ANI_species = []
CSF_species = []
GTDB_species = []
GTDB_genres = []
POP_species = []
sup_genres = []
sup_species = []

# Extract information from each file (the isolate id and clade is contained on its name)
previous_clade = ''
unk_csf_num = len(sources)

for name in all_genomes:

    name = name.replace('.fa', '')
    print(f'Processing {name} ...')

    #__Isolate ID, clade, RefSeq ID, Completeness and Contamination__
    isolate_id = name.split('_')[0]
    isolate_ids.append(isolate_id)
    clade = name.split('_')[1]

    if clade != previous_clade:
        unk_csf_num += 1
        
    clades.append(clade)

    acc = gen_df.loc[gen_df["SAG or Isolate ID"] == isolate_id, "RefSeq Assembly (*IMG Genome ID)"].iloc[0]
    comp = gen_df.loc[gen_df["SAG or Isolate ID"] == isolate_id, "Completeness"].iloc[0]
    cont = gen_df.loc[gen_df["SAG or Isolate ID"] == isolate_id, "Contamination"].iloc[0]

    accs.append(acc)
    comps.append(comp)
    conts.append(cont)

    # __ANI species__
    unk_species_num = len(components) + 1

    A = False
    for i, comp in enumerate(components):
        if name in comp:
            ANI_species.append(i)
            A = True
            break

    if not A:  # Genome does not belong to any specie
        ANI_species.append(unk_species_num)
        unk_species_num += 1

    # __GTDB species and genus__
    GTDB_specie = GTDB_data[name]['s']
    if GTDB_specie == 's__':
        GTDB_species.append('Unk')
    else:
        GTDB_species.append(GTDB_specie.replace('s__', ''))

    GTDB_genus = GTDB_data[name]['g']
    GTDB_genres.append(GTDB_genus.replace('g__', ''))

    # __PopCOGenT species__
    POP_specie = POP_data.loc[POP_data["Strain"] == name, "Main_cluster"].iloc[0]
    POP_species.append(POP_specie)

    # __ConSpeciFix specie__
    unk_csf_num = len(sources) + 1
    csf_clade = CSF_data[clade]

        # Find source in which all the test genomes of the clade belonged to the same specie
    positive_source = [k for k, v in csf_clade.items() if all(x == 1 for x in v)] 
    
    if positive_source:
        CSF_species.append(positive_source)
    else:
        CSF_species.append(unk_csf_num)

    # __Philogenetic information from supplementary data__
    sup_genre = sup_df.loc[sup_df["Subclade Classification"] == clade, "Genus"].iloc[1]
    sup_specie = sup_df.loc[sup_df["Subclade Classification"] == clade, "Species name"].iloc[1]

    sup_genres.append(sup_genre)
    sup_species.append(sup_specie)

    # Remember previous clade
    old_clade = clade

### Final table 
data = {
    "Accession ID": accs,
    "Isolate ID": isolate_ids,
    "Clade": clades, 
    "Completeness": comps,
    "Contamination": conts, 
    "ANI specie": ANI_species,
    "ConSpeciFix specie": CSF_species,
    "PopCOGenT spceie": POP_species,
    "GTDB genus": GTDB_genres,
    "GTDB specie": GTDB_species,
    "Proposed genus": sup_genres,
    "Proposed specie": sup_species
}

df = pd.DataFrame(data)
df.to_csv("genomes_classification.tsv", sep="\t", index=False)
