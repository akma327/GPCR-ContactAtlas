# Author: Anthony Kai Kwang Ma
# Date: 01/05/17
# gen_flareplot_dataset_batch.py

import sys
import csv
import numpy as np

USAGE_STR = """
# Purpose
# Read input flareplot dataset trajectories csv file and output
# bash .sh script that calls simulation_evobundle.py to generate 
# json files for evobundle viewer across full simulations. 

# Usage 
python gen_flareplot_dataset_batch.py <INPUT_TRAJ_TABLE> <FLAREPLOT_DIR> <OUTPUT_FILE> <INTERACTION_TYPE>

# Arguments
# <INPUT_TRAJ_TABLE> Absolute path to the csv file for input trajectories
# <FLAREPLOT_DIR> Path to the directory storing the simulation evobundle jsons
# <OUTPUT_DIR> Absolute path to bash scripts
# <INTERACTION_TYPE> (ie -hbbb) flag to denote specifically which interaction types to 
# display. If no flag provided then defaults to the ['wb', 'wb2', 'lwb', 'lwb2'] collection. 

# Example
INPUT_TRAJ_TABLE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/flareplot_dataset_trajectories.csv"
FLAREPLOT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts"
OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts/batch/gen-flareplot-dataset"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts/batch/gen-flareplot-dataset
python gen_flareplot_dataset_batch.py $INPUT_TRAJ_TABLE $FLAREPLOT_DIR $OUTPUT_DIR -sb

"""

K_MIN_ARG = 4

def write_header(f, INTERACTION_TYPE):
	f.write("#!/bin/bash\n")
	f.write("#Sbatch for: Generate evobundle-simulation flareplots\n")
	f.write("#################\n")
	f.write("#set a job name\n")
	if(INTERACTION_TYPE == ""):
		f.write("#SBATCH --job-name=gen_flareplot_dataset_batch\n")
	else:
		f.write("#SBATCH --job-name=gen_flareplot_dataset_batch_" + INTERACTION_TYPE + "\n")
	f.write("#SBATCH --mail-user=akma327@stanford.edu --mail-type=ALL\n")
	f.write("#################\n")
	f.write("# A file for job output, you can check job progress\n")
	if(INTERACTION_TYPE == ""):
		f.write("#SBATCH --output=gen_flareplot_dataset_batch_PROGRESS.out\n")
	else:
		f.write("#SBATCH --output=gen_flareplot_dataset_batch_" + INTERACTION_TYPE + "_PROGRESS.out\n")
	f.write("#################\n")
	f.write("# A file for errors from the job\n")
	if(INTERACTION_TYPE == ""):
		f.write("#SBATCH --error=gen_flareplot_dataset_batch_ERROR.out\n")
	else:
		f.write("#SBATCH --error=gen_flareplot_dataset_batch_" + INTERACTION_TYPE + "_ERROR.out\n")
	f.write("#################\n")
	f.write("#SBATCH --time=1:00:00\n")
	f.write("#################\n")
	f.write("#SBATCH --partition=rondror\n")
	f.write("#SBATCH --qos=rondror\n")
	f.write("#################\n")
	f.write("#number of nodes you are requesting\n")
	f.write("#SBATCH --tasks=4\n")
	f.write("#SBATCH --ntasks-per-socket=6\n")
	f.write("#SBATCH --mem=10000\n")
	f.write("#################\n\n\n")

def gen_batch_script(INPUT_TRAJ_TABLE, FLAREPLOT_DIR, OUTPUT_DIR, INTERACTION_TYPE=""):
	with open(INPUT_TRAJ_TABLE, 'rU') as f:
		reader = csv.reader(f, delimiter=",")
		d = list(reader)
	darr = np.array(d)
	column_labels = list(darr[0])
	pdb_index = column_labels.index("pdb")
	traj_dir_index = column_labels.index("3.5A_inter_dict_dir")
	output_name_index = column_labels.index("short_exp_cond_description")

	### Extract the pdb, trajectory_directory, output json arguments for each row
	pdb_list = darr[:,pdb_index][1:]
	traj_dir_list = darr[:,traj_dir_index][1:]


	flareplot_subdir = "water-network"
	if(INTERACTION_TYPE != ""):
		flareplot_subdir = INTERACTION_TYPE

	output_file_list = [FLAREPLOT_DIR + "/" + flareplot_subdir + "/" + fp + ".json" for fp in darr[:, output_name_index][1:]]

	### Write out bash script 
	OUTPUT_FILE = OUTPUT_DIR
	if(INTERACTION_TYPE == ""):
		OUTPUT_FILE += "/gen_flareplot_dataset_batch.sbatch"
	else:
		OUTPUT_FILE += "/gen_flareplot_dataset_batch_" + INTERACTION_TYPE + ".sbatch" 
	f = open(OUTPUT_FILE, 'w')
	write_header(f, INTERACTION_TYPE)

	for i, pdb in enumerate(pdb_list):
		traj_dir = traj_dir_list[i]
		output_file = output_file_list[i]
		f.write("EXP_COND_DICT_FOLDER=\"" + traj_dir + "\"\n")
		f.write("EXP_COND_PDB=\"" + pdb + "\"\n")
		f.write("OUTPUT_FILE=\"" + output_file + "\"\n")
		f.write("cd /scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts\n")
		if(INTERACTION_TYPE == ""):
			f.write("python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb &\n\n")
		else:
			f.write("python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb -" + INTERACTION_TYPE + " &\n\n")

	f.write("\n\nwait\n")

if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)
	(INPUT_TRAJ_TABLE, FLAREPLOT_DIR, OUTPUT_DIR) = (sys.argv[1], sys.argv[2], sys.argv[3])
	INTERACTION_TYPE = ""
	if(len(sys.argv) == K_MIN_ARG + 1):
		INTERACTION_TYPE = sys.argv[4].strip("-")
	gen_batch_script(INPUT_TRAJ_TABLE, FLAREPLOT_DIR, OUTPUT_DIR, INTERACTION_TYPE)

