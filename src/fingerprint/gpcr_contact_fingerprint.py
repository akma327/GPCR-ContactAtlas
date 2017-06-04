# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 05/14/17
# gpcr_contact_fingerprint.py

import os
import sys
import json
import csv
import numpy as np
import re
from utils import *

USAGE_STR = """

# Purpose
# Input: GPCR-Contacts for a set of receptors.
# Output: Fingerprint Flareplot format json for visualization. 

# Usage 
# python gpcr_contact_fingerprint.py <OUTPUT_DIR> <INPUT_FILE_1> <INPUT_FILE_2> ...

# Arguments
# <OUTPUT_DIR> Output directory for <interaction_type>.json 
# <INPUT_DIR_1> Input contacts table.txt file for a particular receptor 
# ...
# <INPUT_DIR_N>

# Example
OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/melatonin"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/MTR1A_HUMAN/MTR1A_HUMAN_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/MTR1B_HUMAN/MTR1B_HUMAN_table.txt"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint
python gpcr_contact_fingerprint.py $OUTPUT_DIR $INPUT_FILE_1 $INPUT_FILE_2

OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/drd"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/DRD1_HUMAN/DRD1_HUMAN_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/DRD2_HUMAN/DRD2_HUMAN_table.txt"
INPUT_FILE_3="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/DRD4_HUMAN/DRD4_HUMAN_table.txt"
INPUT_FILE_4="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/homology-models/DRD5_HUMAN/DRD5_HUMAN_table.txt"
INPUT_FILE_5="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/DRD3_HUMAN_3PBL_A/DRD3_HUMAN_3PBL_A_table.txt"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint
python gpcr_contact_fingerprint.py $OUTPUT_DIR $INPUT_FILE_1 $INPUT_FILE_2 $INPUT_FILE_3 $INPUT_FILE_4 $INPUT_FILE_5

OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/opioid"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/OPRM_MOUSE_4DKL_A/OPRM_MOUSE_4DKL_A_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/OPRD_HUMAN_4N6H_A/OPRD_HUMAN_4N6H_A_table.txt"
INPUT_FILE_3="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/OPRK_HUMAN_4DJH_A/OPRK_HUMAN_4DJH_A_table.txt"
INPUT_FILE_4="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/OPRX_HUMAN_4EA3_A/OPRX_HUMAN_4EA3_A_table.txt"
cd /scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint
python gpcr_contact_fingerprint.py $OUTPUT_DIR $INPUT_FILE_1 $INPUT_FILE_2 $INPUT_FILE_3 $INPUT_FILE_4


OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/muscarinic"
INPUT_FILE_1="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/ACM2_HUMAN_3UON_A/ACM2_HUMAN_3UON_A_table.txt"
INPUT_FILE_2="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/ACM3_RAT_4DAJ_A/ACM3_RAT_4DAJ_A_table.txt"
INPUT_FILE_3="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/ACM1_HUMAN_5CXV_A/ACM1_HUMAN_5CXV_A_table.txt"
INPUT_FILE_4="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs/ACM4_HUMAN_5DSG_A/ACM4_HUMAN_5DSG_A_table.txt"
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
	frameDict = {}
	edge_dict = {}
	for i in range(len(sorted_input_files)):
		p, label = sorted_input_files[i], labels[i]
		frameDict[i] = label 

		f = open(p, 'r')
		for line in f:
			gpcrdb1, gpcrdb2, res1, res2, itype = line.strip().split("\t")
			if(itype not in edge_dict):
				edge_dict[itype] = {}
			if((gpcrdb1, gpcrdb2) not in edge_dict[itype]):
				edge_dict[itype][(gpcrdb1, gpcrdb2)] = set()
			edge_dict[itype][(gpcrdb1, gpcrdb2)].append(i)

	return edge_dict, frameDict



def fingerprint(OUTPUT_DIR, INPUT_FILES):
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)

	sorted_input_files, labels = sort_paths(INPUT_FILES)
	edge_dict, frameDict = get_edge_and_frameDict(sorted_input_files, labels)
	for itype in edge_dict:
		out_path = OUTPUT_DIR + "/" + itype + ".json"
		edges = []
		for edge in edge_dict[itype]:
			name1, name2 = edge
			frames = sorted(list(edge_dict[itype][edge]))
			edges.append({"name1": name1, "name2": name2, "frames": frames})
		
		json_dict = partial_input
		json_dict["edges"] = edges
		json_dict["frameDict"] = frameDict
		with open(out_path, 'w') as outfile:
			json.dump(json_dict, outfile)




if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)
	OUTPUT_DIR = sys.argv[1]
	INPUT_FILES = sys.argv[2:]
	fingerprint(OUTPUT_DIR, INPUT_FILES)
