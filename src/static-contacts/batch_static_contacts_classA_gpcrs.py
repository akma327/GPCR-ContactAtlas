# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 04/29/17
# batch_static_contacts_classA_gpcrs.py

import os
import sys
import glob 
import json
from utils import *

USAGE_STR = """
# Purpose
# Generate .sh script to call StaticInteractionCalculator.py upon every available
# homology model or existing crystal structure. 

# Usage 
# python batch_static_contacts_classA_gpcrs.py 

"""


HADDED_CLASS_A_GPCR_PATH="/scratch/PI/rondror/akma327/GPCRContacts/data/classA-gpcr-pdbs"
STATIC_CONTACTS_CLASS_A_GPCR_PATH="/scratch/PI/rondror/akma327/GPCRContacts/data/static-contacts/classA-gpcr-pdbs"
PDB_TO_LIGAND_PATH = "/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/ligand-wetness/watermarks/tables/classA_GPCRs_known_ligands.tsv"

def pdb_to_ligand():
	"""
		Create mapping from pdb to ligand 
	"""
	pdb_to_ligand_dict = {}
	f = open(PDB_TO_LIGAND_PATH, 'r')
	header = f.readline()
	for line in f:
		uniprot, class_code, pdb, ligand = line.strip().split("\t")
		pdb_to_ligand_dict[pdb] = ligand

	return pdb_to_ligand_dict

def compute_classA_gpcr_contacts():
	"""
		Compute static contacts for all class A gpcr crystals
	"""
	pdb_to_ligand_dict = pdb_to_ligand()
	hadded_gpcr_paths = glob.glob(HADDED_CLASS_A_GPCR_PATH + "/*")
	for i, gpcr_pdb_path in enumerate(hadded_gpcr_paths):
		# if(i > 0): break
		print(gpcr_pdb_path)
		pdb, chain = gpcr_pdb_path.split("/")[-1].split("_")[0:2]
		uniprot = getUniprotCode(pdb)
		identifier = uniprot + "_" + pdb + "_" + chain

		output_dir = STATIC_CONTACTS_CLASS_A_GPCR_PATH + "/" + identifier
		output_contacts = output_dir + "/" + identifier + "_contacts.txt"

		ligand = pdb_to_ligand_dict[pdb]
		calc_command = "python StaticInteractionCalculator.py " + gpcr_pdb_path + " " + output_contacts + " -interlist -sb -pc -ps -ts -vdw -hbbb -hbsb -hbss -wb -wb2 -lwb -lwb2 -hlb -hls -ligand " + ligand
		os.chdir("/scratch/PI/rondror/akma327/DynamicNetworks/src/interaction-geometry")
		os.system(calc_command)


		### Convert the raw contacts.txt to UNIPROT_table.txt file 
		RESI_TO_GPCRDB = genGpcrdbDict(uniprot)
		output_table = output_dir + "/" + identifier + "_table.txt"
		f = open(output_contacts, 'r')
		fw = open(output_table, 'w')

		contact_type_to_edges = {}
		for line in f:
			atom_pair, contact_type = line.strip().split("@-")
			atom1, atom2 = atom_pair.split(" -- ")
			resi1, resi2 = atom1.split("-")[0], atom2.split("-")[0]
			gpcrdb1, gpcrdb2 = getGPCRDB(resi1, RESI_TO_GPCRDB), getGPCRDB(resi2, RESI_TO_GPCRDB)
			if(gpcrdb1 == "None" or gpcrdb2 == "None"): continue 
			if(not ligOrInTM(gpcrdb1) or not ligOrInTM(gpcrdb2)): continue

			gpcrdb1, gpcrdb2, resi1, resi2 = orderGpcrdbsAndResi(gpcrdb1, gpcrdb2, resi1, resi2)
			key = (gpcrdb1, gpcrdb2, resi1, resi2)
			if(contact_type not in contact_type_to_edges):
				contact_type_to_edges[contact_type] = set()
				contact_type_to_edges[contact_type].add(key)
			else:
				if(key not in contact_type_to_edges[contact_type]):
					contact_type_to_edges[contact_type].add(key)

		for contact_type in contact_type_to_edges:
			for gpcrdb1, gpcrdb2, resi1, resi2 in contact_type_to_edges[contact_type]:
				fw.write(gpcrdb1 + "\t" + gpcrdb2 + "\t" + resi1 + "\t" + resi2 + "\t" + contact_type + "\n")

		fw.close()

		### Convert UNIPROT_table.txt to corresponding flareplot json files 
		f2 = open(output_table, 'r')
		contact_type_to_edges = {}
		for line in f2:
			gpcrdb1, gpcrdb2, resi1, resi2, contact_type = line.strip().split("\t")
			gpcrdb1, gpcrdb2 = orderGpcrdbs(gpcrdb1, gpcrdb2)
			key = (gpcrdb1, gpcrdb2)
			if(contact_type not in contact_type_to_edges):
				contact_type_to_edges[contact_type] = set()
				contact_type_to_edges[contact_type].add(key)
			else:
				if(key not in contact_type_to_edges[contact_type]):
					contact_type_to_edges[contact_type].add(key)

		### Generate a json for every interaction type
		for interaction_type in contact_type_to_edges:
			output_json_path = output_dir + "/" + identifier + "_" + interaction_type + ".json"
			json_dict = node_dict
			edges = []
			for name1, name2 in contact_type_to_edges[interaction_type]:
				edges.append({"name1": name1, "name2": name2, "frames": [0]})
			json_dict["edges"] = edges

			with open(output_json_path, 'w') as f:
				json.dump(json_dict, f)


if __name__ == "__main__":
	compute_classA_gpcr_contacts()



