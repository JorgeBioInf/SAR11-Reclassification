"""
CSF_clades_analysis.py
----------------------
Runs ConSpeciFix analysis between all test groups and source groups

Author: Jorge Marcos Fernández
Date: 2025-12-18
Version: 1.0

Usage:
    python CSF_clades_analysis.py genomes_list.txt source_dir1 source_dir2 ... source_dirn

Output:
    - CSF_clades_results.json dictionary with analysis results

Dependencies:
    - conda
    - conda environment "CSF" with python 2.7.
    - ConSpeciFix
    - sys
    - json
    - os
    - shutil
    - subprocess
    - collections
    
Notes:
    - Requires a .txt file with a list containing the names of all test genomes 
    - Requires test groups folders (divided by clades) with names clade_{num}/ in current directory
    - May require change ConSpeciFix runner_personal.py file location in the code
    - Recommended to run in background
"""


# -- PACKAGES --
import subprocess
import sys
import os
import shutil
import random
import json
from collections import defaultdict

# -- ARGUMENTS CHECK --
all_genomes_file = sys.argv[1]

try:
    with open(all_genomes_file, 'r') as file:
        all_genomes = [line.strip() for line in file]
except Exception as e:
    print(f'Error reading file with all genomes: {e}')
    sys.exit(1)

srcs = [s for s in sys.argv[2:len(sys.argv)]]

for src in srcs:
    if not os.path.exists(src):
        print(f'Error: no folder named {src}')
        print('Please, provide a valid input')
        sys.exit(1)

# -- FUNCTIONS --
def run_conspecific(subfolder_path):
    """
    Executes ConSpeciFix Python 2.7 using conda environment 'CSF'
    """

    command = [
        "conda", "run", "-n", "CSF",
        "python",
        "/home/estudiante2/JMF/ConSpeciFix/ConSpeciFix-1.3.0/database/runner_personal.py", # Change to ConSpeciFix location
        subfolder_path
    ]

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )

    except subprocess.CalledProcessError as e:
        print("ERROR running ConSpeciFix")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)


def parse_results(results_dir, genome):
    """Parses a results.txt file from ConSpeciFix and checks whether test genome belong to same species or not"""
    # Read results file
    with open(results_dir, 'r') as file:
        text = file.readlines()

    # Parse text and store same species and different species
    species = False
    not_species = False

    same_species = []
    diff_species = []

    for line in text:
        line = line.rstrip('\n')

        if line == 'The following strains are members of the species:':
            species = True
            continue

        if line == 'The following strains were determined to NOT be a member of the species:':
            species = False
            continue

        if species and line:
            same_species.append(line)

    print(f'\n\nProcessing test genome {genome} ...')

    # Check if same species
    if genome in same_species:
        return True
    else:
        return False
    
def extract_and_copy_gno2(new_dir, out_dir, new_name):
    "Searches for gno2 plot and stores it in directory results_plots with a new name"

    gno2_path = os.path.join(new_dir, "_conspecifix/database/User_spec/gno2.png")

    dst_path = os.path.join(out_dir, f"{new_name}")

    try:
        shutil.copy(gno2_path, dst_path)
        print(f"Saved plot: {dst_path}")
    except Exception as e:
        print(f'Error copying gon2.png plot: {e}')


# -- MAIN PROGRAM --

# Get all possible clades
all_clades = []
for genome in all_genomes:
    all_clades.append(genome.split('_')[1].replace('.fa', ''))

unique_clades = set(all_clades)

# Create directory to store interesting plots
out_dir = "CSF_results_and_plots"
os.makedirs(out_dir, exist_ok=True)

# Structure with final results
all_output = {}

for clade in unique_clades:
    folder_name = f'clade_{clade}'
    test_num = 2

    clade_dict = defaultdict(list)

    # Check if groups folder exists
    if not os.path.isdir(folder_name): 
        print(f'No folder for clade {clade}')
        continue

    # Select random genome from the folder name
    files = [f for f in os.listdir(folder_name) if os.path.isfile(os.path.join(folder_name, f))] 

    if len(files) < test_num:
        test_num = len(files)

    for genome in random.sample(files, k=test_num):   

        genome_name = genome.replace('.fa', '').replace('anotated_', '')
        genome_path = os.path.join(folder_name, genome)

        # Test with each source 
        for source_num, source_folder in enumerate(srcs):
            source_num_str = str(source_num + 1)
            analysis_folder_name = f'analysis_folder_{clade}_{genome_name}_s{source_num + 1}'

            os.makedirs(analysis_folder_name)
            analysis_folder_abspath = os.path.abspath(analysis_folder_name)

            # Copy test genome
            shutil.copy(genome_path, os.path.join(analysis_folder_name, genome))

            # Copy source genomes
            for f in os.listdir(source_folder):
                src_f = os.path.join(source_folder, f)
                dst_f = os.path.join(analysis_folder_name, f)
                if os.path.isfile(src_f):
                    shutil.copy(src_f, dst_f)

            # Run ConSpeciFix
            print(f'\nRunning ConSpeciFix between {genome_name} in clade {clade} and source {source_num + 1}')
            run_conspecific(analysis_folder_abspath)

            # Parse results.txt file 
            results_path = os.path.join(analysis_folder_name, "results.txt")

            if not os.path.isfile(results_path):
                print("Error: no file results.txt retrieved for the analysis\n")
                shutil.rmtree(analysis_folder_name)
                continue  

            # Extract results.txt gno2.png plot and store it in output folder
            analysis = f'{genome}_{clade}_s{source_num}' 
            resultsname = f'{analysis}_results.txt'
            results_out_path = os.path.join(out_dir, resultsname)
            shutil.copy(results_path, results_out_path)

            plotname = f'{analysis}_gno2.png'
            extract_and_copy_gno2(analysis_folder_name, out_dir, plotname)

            # Search if test genome belongs to same species as source
            species = parse_results(results_path, genome)
            if species:
                print(f'{genome_name} belongs to same specie as group {source_folder}!')
                clade_dict[source_num_str].append(1)
            else:
                print(f'{genome_name} DOES NOT belong to same specie as group {source_folder}!!')
                clade_dict[source_num_str].append(0)

            shutil.rmtree(analysis_folder_name)

    all_output[clade] = clade_dict

    print('')
    print('')

# Store final results
with open('CSF_clades_results.json', 'w') as file:
    json.dump(all_output, file, indent = 4)

print('')
print('')
print('\nAnalysis finished!')
