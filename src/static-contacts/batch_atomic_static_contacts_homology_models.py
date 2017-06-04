# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 04/29/17
# batch_atomic_static_contacts_homology_models.py

import os
import sys
import glob 
import json
from utils import *

USAGE_STR = """
# Purpose
# Generate .sh script to call StaticInteractionCalculator.py upon every available
# homology model or existing crystal structure. Includes atom identity 

# Usage 
# python batch_atomic_static_contacts_homology_models.py 

"""

HADDED_HOMOLOGY_MODEL_PATH="/scratch/PI/rondror/akma327/GPCRContacts/data/structures/hadded-homology-models"
STATIC_CONTACTS_HM_PATH="/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/homology-models"


def compute_homology_model_contacts():
	"""
		Compute static contacts for all homology models
	"""
	print("compute_homology_model_contacts ...")

	hadded_hm_paths = glob.glob(HADDED_HOMOLOGY_MODEL_PATH + "/*")
	for i, hm_pdb_path in enumerate(hadded_hm_paths):
		# if(i > 1): break
		print(hm_pdb_path)
		uniprot = hm_pdb_path.split("/")[-1].split(".")[0]
		output_dir = STATIC_CONTACTS_HM_PATH + "/" + uniprot
		output_contacts = output_dir + "/" + uniprot + "_contacts.txt"
		calc_command = "python StaticInteractionCalculator_water_indexed_crys.py " + hm_pdb_path + " " +  output_contacts + " -interlist -all"
		os.chdir("/scratch/PI/rondror/akma327/DynamicNetworks/src/interaction-geometry")
		os.system(calc_command)


		### Convert the raw contacts.txt to UNIPROT_table.txt file 
		RESI_TO_GPCRDB = genGpcrdbDict(uniprot)
		output_table = output_dir + "/" + uniprot + "_table.txt"
		f = open(output_contacts, 'r')
		fw = open(output_table, 'w')

		contact_type_to_edges = {}
		for line in f:
			atoms, contact_type = line.strip().split("@-")
			atoms = atoms.split(" -- ")

			### Clean off chain data 
			atoms = [a.split("_")[0] for a in atoms]
			print(atoms)

			resatom1, resatom2 = atoms[0].split("-"), atoms[1].split("-")
			if(len(resatom1) !=2 or len(resatom2) != 2): continue

			resi1, atom1 = resatom1
			resi2, atom2 = resatom2

			gpcrdb1, gpcrdb2 = getGPCRDB(resi1, RESI_TO_GPCRDB), getGPCRDB(resi2, RESI_TO_GPCRDB)
			if(gpcrdb1 == "None" or gpcrdb2 == "None"): continue 
			if(not ligOrInTM(gpcrdb1) or not ligOrInTM(gpcrdb2)): continue

			if(flipGpcrdbs(gpcrdb1, gpcrdb2) == True):
				if(len(atoms) == 2):
					atoms = [atoms[1], atoms[0]]
				elif(len(atoms) == 3):
					atoms = [atoms[1], atoms[0], atoms[2]]
				elif(len(atoms) == 4):
					atoms = [atoms[1], atoms[0], atoms[3], atoms[2]]
				key = (gpcrdb2, gpcrdb1, ":".join(atoms))
			else:
				key = (gpcrdb1, gpcrdb2, ":".join(atoms))

			if(contact_type not in contact_type_to_edges):
				contact_type_to_edges[contact_type] = set()
				contact_type_to_edges[contact_type].add(key)
			else:
				if(key not in contact_type_to_edges[contact_type]):
					contact_type_to_edges[contact_type].add(key)


		for contact_type in contact_type_to_edges:
			for gpcrdb1, gpcrdb2, atom_list in contact_type_to_edges[contact_type]:
				fw.write(gpcrdb1 + "\t" + gpcrdb2 + "\t" + atom_list + "\t" + contact_type + "\n")

		fw.close()
		

		### Convert UNIPROT_table.txt to corresponding flareplot json files 
		f2 = open(output_table, 'r')
		contact_type_to_edges = {}
		for line in f2:
			gpcrdb1, gpcrdb2, atom_list, contact_type = line.strip().split("\t")
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
			output_json_path = output_dir + "/" + uniprot + "_" + interaction_type + ".json"
			json_dict = node_dict
			edges = []
			for name1, name2 in contact_type_to_edges[interaction_type]:
				edges.append({"name1": name1, "name2": name2, "frames": [0]})
			json_dict["edges"] = edges

			with open(output_json_path, 'w') as f:
				json.dump(json_dict, f)


if __name__ == "__main__":
	compute_homology_model_contacts()


