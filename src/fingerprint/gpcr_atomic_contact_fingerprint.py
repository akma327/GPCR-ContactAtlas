# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 05/14/17
# gpcr_atomic_contact_fingerprint.py

import os
import sys
import json
import csv
import numpy as np
import re
from utils import *

USAGE_STR = """

# Purpose
# Input: GPCR-Contacts at atomic resolution for a set of receptors.
# Output: Fingerprint Flareplot format json for visualization. 

# Usage 
# python gpcr_contact_fingerprint.py <OUTPUT_DIR> <INPUT_FILE_1> <INPUT_FILE_2> ...

# Arguments
# <OUTPUT_DIR> Output directory for <interaction_type>.json 
# <INPUT_DIR_1> Input contacts table.txt file for a particular receptor 
# ...
# <INPUT_DIR_N>

# Example
OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/opioid"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/OPRM_MOUSE_4DKL_A/OPRM_MOUSE_4DKL_A_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/OPRD_HUMAN_4N6H_A/OPRD_HUMAN_4N6H_A_table.txt"
INPUT_FILE_3="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/OPRK_HUMAN_4DJH_A/OPRK_HUMAN_4DJH_A_table.txt"
INPUT_FILE_4="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/OPRX_HUMAN_4EA3_A/OPRX_HUMAN_4EA3_A_table.txt"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint
python gpcr_contact_fingerprint.py $OUTPUT_DIR $INPUT_FILE_1 $INPUT_FILE_2 $INPUT_FILE_3 $INPUT_FILE_4


OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/muscarinic"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/ACM2_HUMAN_3UON_A/ACM2_HUMAN_3UON_A_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/ACM3_RAT_4DAJ_A/ACM3_RAT_4DAJ_A_table.txt"
INPUT_FILE_3="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/ACM1_HUMAN_5CXV_A/ACM1_HUMAN_5CXV_A_table.txt"
INPUT_FILE_4="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs/ACM4_HUMAN_5DSG_A/ACM4_HUMAN_5DSG_A_table.txt"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint
python gpcr_contact_fingerprint.py $OUTPUT_DIR $INPUT_FILE_1 $INPUT_FILE_2 $INPUT_FILE_3 $INPUT_FILE_4

"""

K_MIN_ARG = 3

def sort_paths(INPUT_FILES):
	path_comp = [p.split("/") for p in INPUT_FILES]
	sorted_input_files = ["/".join(sp) for sp in sorted(path_comp, key=lambda x: x[-1])]
	labels = [sp[-1].strip("_table.txt") for sp in sorted(path_comp, key = lambda x: x[-1])]
	return sorted_input_files, labels


def get_edge_and_frameDict(sorted_input_files, labels):
	"""
		sorted_input_files: Tables containing the gpcrdb pairs and atomic resolution interaction data 
		labels: GPCR Receptor name 
	"""
	frameDict = {}
	edge_dict = {}

	# residue_to_atomic_contacts = {itype: {gpcrdb1:gpcrdb2: {receptor_index1: [list of atomic resolution contacts], ...}}}
	# {"hbbb": {3x53:3x54: {0:[ALA168-N:VAL169-N, ALA168-O:VAL169-N,...], 1: []}}}
	residue_to_atomic_contacts = {} 
	for i in range(len(sorted_input_files)):
		p, label = sorted_input_files[i], labels[i]
		frameDict[i] = label 

		f = open(p, 'r')
		for line in f:
			gpcrdb1, gpcrdb2, atom_list, itype = line.strip().split("\t")
			if(itype not in edge_dict):
				edge_dict[itype] = {}
			if((gpcrdb1, gpcrdb2) not in edge_dict[itype]):
				edge_dict[itype][(gpcrdb1, gpcrdb2)] = []
			edge_dict[itype][(gpcrdb1, gpcrdb2)].append(i)


			### Atomic level contacts encoding into flareplot json
			gpcrdb_pair = (gpcrdb1, gpcrdb2)
			if(itype not in residue_to_atomic_contacts):
				residue_to_atomic_contacts[itype] = {}
			if(gpcrdb_pair not in residue_to_atomic_contacts[itype]):
				residue_to_atomic_contacts[itype][gpcrdb_pair] = {}
			if(atom_list not in residue_to_atomic_contacts[itype][gpcrdb_pair]):
				residue_to_atomic_contacts[itype][gpcrdb_pair][i] = []
			residue_to_atomic_contacts[itype][gpcrdb_pair][i].append(atom_list)

	return edge_dict, frameDict, residue_to_atomic_contacts



def fingerprint(OUTPUT_DIR, INPUT_FILES):
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)

	sorted_input_files, labels = sort_paths(INPUT_FILES)
	edge_dict, frameDict, residue_to_atomic_contacts = get_edge_and_frameDict(sorted_input_files, labels)
	for itype in edge_dict:
		out_path = OUTPUT_DIR + "/" + itype + ".json"
		edges = []
		for edge in edge_dict[itype]:
			name1, name2 = edge
			frames = edge_dict[itype][edge]
			edges.append({"name1": name1, "name2": name2, "frames": frames})
		
		json_dict = partial_input
		json_dict["edges"] = edges
		json_dict["frameDict"] = frameDict

		### Atomic level contacts encoding into flareplot json
		json_dict["residue_to_atomic_contacts"] = residue_to_atomic_contacts[itype]


		with open(out_path, 'w') as outfile:
			json.dump(json_dict, outfile)




if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)
	OUTPUT_DIR = sys.argv[1]
	INPUT_FILES = sys.argv[2:]
	fingerprint(OUTPUT_DIR, INPUT_FILES)
