# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 06/06/17
# simulation_evobundle.py 

import os
import sys 
import json 
import glob 
import ast 
import pickle
import time
from utils import * 

USAGE_STR = """

# Purpose 
# Given a path to the directory containing interaction dictionaries for a particular experiment
# condition, generate the evobundle json file to play all these simulate these interactions
# on a bundle plot. User has an option to view nodes as GPCRDB numbers or original amino acids. 

# Usage 
# python simulation_evobundle.py <EXP_COND_DICT_FOLDER> <EXP_COND_PDB> <OUTPUT_FILE> <NODE_FLAG> <INTERACTION_TYPE>

# Arguments
# <EXP_COND_DICT_FOLDER> Absolute path to the experiment replicate dictionary folder. 
# <EXP_COND_PDB> PDB Code for the experiment condition to generate simulation evobundle json for. 
# <OUTPUT_FILE> Absolute path to output simulation evobundle json file. 
# <NODE_FLAG> 
#	gpcrdb To show the gpcrdb numbers (ie 5x50)
#	residue To show explicit resname and resnum (ie ASP130)
# <INTERACTION_TYPE> Optional flag to denote a single interaction type to compute for. Defaults 
# to ['wb', 'wb2', 'lwb', 'lwb2'] collection if no flag provided. 

# Example 
EXP_COND_DICT_FOLDER="/scratch/PI/rondror/akma327/DynamicNetworks/data/DynamicNetworksOutput/InteractionOutputDictionary/MOR-active-noNb-BU72-nature2015/condition-nb-interface-restrained-BU72/rep_1/050516_NONHB_DEFAULT_HB_3.5_70"
EXP_COND_PDB="5C1M"
OUTPUT_FILE="/scratch/PI/rondror/akma327/GPCRContacts/data/dynamic-contacts/MOR:active-noNb:BU72-nB_interface_restrained:3.0us:312.18K.json"
cd /scratch/PI/rondror/akma327/DynamicNetworks/src/analysis/interaction-footprint-analysis
python simulation_evobundle.py $EXP_COND_DICT_FOLDER $EXP_COND_PDB $OUTPUT_FILE gpcrdb

"""

K_MIN_ARG = 5

# INTERACTION_FILTER = ['sb', 'pc', 'ps', 'ts', 'vdw', 'hbbb', 'hbsb', 'hbss', 'wb', 'wb2']
# INTERACTION_FILTER = ["sb", 'pc', 'ps', 'ts']
INTERACTION_FILTER = ["wb", "wb2", 'lwb', 'lwb2']
# INTERACTION_FILTER = ['lwb', 'lwb2']

SIMULATION_TABLE="/scratch/PI/rondror/akma327/GPCRContacts/src/dynamic-contacts/simulation_table.tsv"


def residuesTooClose(gpcrdb1, gpcrdb2, contact_type):
	"""
		Don't count inter-residue interactions or van der waals that
		occur between residues within 4 indices apart 
	"""
	if("LIG" in gpcrdb1 or "LIG" in gpcrdb2): return False
	tm1, index1 = map(int, gpcrdb1.split("x"))
	tm2, index2 = map(int, gpcrdb2.split("x"))
	if(tm1 != tm2): return False
	if("vdw" in contact_type):
		if(index1 > index2 - 4 and index1 < index2 + 4): return True
	else:
		if(index1 == index2): return True
	return False

# Filters list of interaction dictionary file paths and only keeps the ones 
# that are part of the list of desired interaction type
def interactionFilter(interaction_dict_files, INTERACTION_FILTER):
	filtered_files = set()
	for dict_file in interaction_dict_files:
		for interaction_flag in INTERACTION_FILTER:
			dict_file_suffix = dict_file.split("/")[-1]
			if(flag_to_dict_file[interaction_flag] == dict_file_suffix):
				filtered_files.add((dict_file,interaction_flag))
	return list(filtered_files)


# Get interaction dictionary files 
def getInteractionDictFiles(EXP_COND_DICT_FOLDER):
	interaction_dict_files = glob.glob(EXP_COND_DICT_FOLDER + "/*dict.txt")
	return interactionFilter(interaction_dict_files, INTERACTION_FILTER)


# Fill json dictionary with list of interaction keys 
def fillInteractionKeys(interaction_dict_files, EXP_COND_PDB_DICT, NODE_FLAG):
	# {(gpcrdb1, gpcrdb2): Longest time point list} for gpcrdb case 
	# {(res1, res2): Longest time point list} for explicit amino acid case 
	interactionDict = {} 
	#### Iterate through all interaction file types 
	for interaction_dict_file, interaction_flag in interaction_dict_files:
		print(interaction_dict_file, interaction_flag)

		f = open(interaction_dict_file, 'r')
		totalFrames = None 
		for line in f: 
			### Retrieve total frame count 
			if("TotalFrames" in line): totalFrames = int(line.split("TotalFrames:")[1].strip())
			### Process single interaction 
			if("--" in line):
				resAtomPair, tpstr = line.split("~")
				res1, res2 = extractResPair(resAtomPair)
				gpcrdb1, gpcrdb2 = extractGpcrdbPair(resAtomPair, EXP_COND_PDB_DICT)
				if(gpcrdb1 == "None" or gpcrdb2 == "None"): continue 
				if(gpcrdb1 == "------" or gpcrdb2 == "------"): continue
				if(NODE_FLAG == "gpcrdb"): key1, key2 = (gpcrdb1, gpcrdb2, interaction_flag), (gpcrdb2, gpcrdb1, interaction_flag)
				if(NODE_FLAG == "residue"): key1, key2 = (res1, res2, interaction_flag), (res2, res1, interaction_flag)
				if(key1 not in interactionDict and key2 not in interactionDict):
					interactionDict[key1] = ast.literal_eval(tpstr)
				else:
					numTimepoints = len(tpstr.split(","))
					if(key1 in interactionDict):
						if(numTimepoints > len(interactionDict[key1])): interactionDict[key1] = ast.literal_eval(tpstr)
					else: 
						if(numTimepoints > len(interactionDict[key2])): interactionDict[key2] = ast.literal_eval(tpstr)
	return totalFrames, interactionDict


def simulation_id_to_ligand():
	"""
		Create mapping from simulation identifier to ligand
	"""
	f = open(SIMULATION_TABLE, 'r')
	simulation_id_to_lig = {}
	for line in f:
		uniprot, pdb, simulation_id, ligand = line.strip().split("\t")
		simulation_id_to_lig[simulation_id] = ligand 

	return simulation_id_to_lig

def simulationEvobundleDriver(EXP_COND_DICT_FOLDER, EXP_COND_PDB, OUTPUT_FILE, NODE_FLAG):

	tic = time.clock()
	OUT_DIR = "/".join(OUTPUT_FILE.split("/")[:-1])
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	print("Compute Simulation EvoBundle Json: " + EXP_COND_DICT_FOLDER)
	EXP_COND_PDB_DICT = genResidueToGpcrdbTable(EXP_COND_PDB)
	interaction_dict_files = getInteractionDictFiles(EXP_COND_DICT_FOLDER)
	totalFrames, interactionDict = fillInteractionKeys(interaction_dict_files, EXP_COND_PDB_DICT, NODE_FLAG)

	### Write to json file 
	simulation_id_to_lig = simulation_id_to_ligand()
	simulation_id = OUTPUT_FILE.split("/")[-1].strip(".json")
	ligand = "LIG"
	if(simulation_id in simulation_id_to_lig):
		ligand = simulation_id_to_lig[simulation_id]

	json_dict = partial_input
	edges = []
	for key in interactionDict:
		gpcrdb1, gpcrdb2, interaction_type = key 
		tp = interactionDict[key]
		name1 = gpcrdb1
		name2 = gpcrdb2

		# No intra-residue interactions and VDW interactions within 4 residues distance
		contact_type = INTERACTION_FILTER[0]
		if(residuesTooClose(gpcrdb1, gpcrdb2, contact_type)): continue
		if(name1 == "LIGxLIG"):
			name1 = "LIGx" + ligand 
		if(name2 == "LIGxLIG"):
			name1 = "LIGx" + ligand



		edges.append({"name1": name1, "name2": name2, "frames": tp})

	json_dict["edges"] = edges 

	with open(OUTPUT_FILE, 'w') as outfile:
		json.dump(json_dict, outfile)

	toc = time.clock()
	print("Time: " + str(toc - tic))




if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(1)
	(EXP_COND_DICT_FOLDER, EXP_COND_PDB, OUTPUT_FILE, NODE_FLAG) = (sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	if(len(sys.argv) == K_MIN_ARG + 1):
		INTERACTION_FILTER = [sys.argv[5].strip("-")]
		
	simulationEvobundleDriver(EXP_COND_DICT_FOLDER, EXP_COND_PDB, OUTPUT_FILE, NODE_FLAG)
