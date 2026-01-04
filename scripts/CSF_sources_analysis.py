"""
CSF_sources_analysis.py
----------------------
Runs ConSpeciFix analysis among source groups

Author: Jorge Marcos Fernández
Date: 2025-12-15
Version: 1.0

Usage:
    python CSF_sources_analysis.py dir1 dir2 ... dirn

Output:
    - CSF_source_results.json dictionary with analysis results

Dependencies:
    - conda
    - conda environment "CSF" with python 2.7.
    - ConSpeciFix
    - sys
    - json
    - os
    - shutil
    - subprocess
    
Notes:
    - Only needed if many source groups are formed during ANI clustering
    - Can process as many directories as given. Note that computing time may considerably increase if a great number is offered
    - May require change ConSpeciFix runner_personal.py file location in the code
    - Recommended to run in background  
"""

# -- PACKAGES --
import subprocess
import sys
import os
import random
import shutil
import json

# -- ARGUMENTS CHECK --
if len(sys.argv) <= 2:
    print('Use: CSF_sources_analysis.py dir1 dir2 ... dirn')
    sys.exit(1)

srcs = [s for s in sys.argv[1:]]

for src in srcs:
    if not os.path.exists(src):
        print(f'Error: no folder named {src}')
        print('Please, provide a valid input')
        sys.exit(1)

# -- FUNCTIONS --
def run_conspecific(subfolder_path, log_file="conspecific_output.txt"):
    """
    Executes ConSpeciFix Python 2.7 using conda environment 'CSF'
    and saves all stdout and stderr to a log file.
    """

    command = [
        "conda", "run", "-n", "CSF",
        "python",
        "/home/estudiante2/JMF/ConSpeciFix/ConSpeciFix-1.3.0/database/runner_personal.py",
        subfolder_path
    ]
    print(subfolder_path)
    with open(log_file, "w") as f:  
        try:
            result = subprocess.run(
                command,
                stdout=f,         
                stderr=subprocess.STDOUT,  
                text=True,
                check=True
            )
            print(f"ConSpeciFix finished. Output saved to {log_file}")

        except subprocess.CalledProcessError as e:
            print(f"ERROR running ConSpeciFix. Check {log_file} for details.")


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
            not_species = True
            continue

        if species and line:
            same_species.append(line)

        if not_species and line:
            diff_species.append(line)

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

# Select a random test genome for each folder
random_candidates = {}  # Dictionary with the relation random_genome (selected from source folder) 
                        # - source_index (number associated to the specific source folder)
for idx, src in enumerate(srcs):
    files = [f for f in os.listdir(src) if os.path.isfile(os.path.join(src, f))] # Store all files of the folder
    random_test = random.choice(files)
    random_candidates[random_test] = idx + 1

print('The following genomes have been selected as random candidates:')
for genome, idx in random_candidates.items():
    print(f'{genome} from source folder {idx}')
print('')


# Create directory to store interesting plots
out_dir = "CSF_results_and_plots"
os.makedirs(out_dir, exist_ok=True)

# Execute ConSpeciFix for each test-source pair
results_dict = {}

for genome, idx in random_candidates.items():
    for element in range(len(srcs)):
        i = element + 1
        if idx != i: # Avoid evaluating each test over its own group
            print(f'\nEvaluating genome from source folder {idx} with group {i}')
            source_folder = srcs[element]  # Source genomes path to compare
            original_folder = srcs[idx - 1] # Test genome's source group path
            genome_path = os.path.join(original_folder, genome) # Genome full path

            # Create new directory
            new_dir = f'CSF_{idx}_{i}'
            os.makedirs(new_dir, exist_ok=True)
            abs_path = os.path.abspath(new_dir)                        
            # Copy genomes of source folder
            for f in os.listdir(source_folder):
                src_f = os.path.join(source_folder, f)
                dst_f = os.path.join(new_dir, f)
                if os.path.isfile(src_f):
                    shutil.copy(src_f, dst_f)

            # Copy test genome
            shutil.copy(genome_path, os.path.join(new_dir, genome))

            # Run ConSpeciFix
            print('Runnning ConSpeciFix ...')
            print('')
            run_conspecific(abs_path)

            # Parse results.txt file and check if test genome belongs to same species
            results_path = os.path.join(abs_path, "results.txt")

            if not os.path.isfile(results_path):
                print("Error: no file results.txt retrieved for the analysis")
                continue

            species = parse_results(results_path, genome)

            # Store results
            term = f'{idx}-{i}'
            if species:
                results_dict[term] = 'YES'
            else:
                results_dict[term] = 'NO'

            # Extract results.txt gno2.png plot and store it in output folder
            resultsname = f'{term}_results.txt'
            results_out_path = os.path.join(out_dir, resultsname)
            shutil.copy(results_path, results_out_path)

            plotname = f'{term}_gno2.png'
            extract_and_copy_gno2(new_dir, out_dir, plotname)

            # Delete temporary folder
            shutil.rmtree(new_dir)

            print('')
            print('')


# Store final results
with open('CSF_source_results.json', 'w') as file:
    json.dump(results_dict, file, indent = 4)
