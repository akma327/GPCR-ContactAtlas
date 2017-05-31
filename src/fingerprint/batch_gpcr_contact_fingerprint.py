# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 05/18/17
# batch_gpcr_contact_fingerprint.py

import os
import sys
import json
import glob


USAGE_STR = """

# Purpose
# Call gpcr_contact_fingerprint.py function on every batch of receptor uniprots

# Usage 
# python gpcr_contact_fingerprint.py <OUTPUT_DIR>

# Arguments
# <OUTPUT_DIR> Absolute path to the output directory that will contain all the flareplot jsons for each contact type
# for every cluster of uniprots 

# Example
OUTPUT_DIR="/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/051817"
python batch_gpcr_contact_fingerprint.py $OUTPUT_DIR

"""

K_MIN_ARG = 2

gpcr_tree_path = "/scratch/PI/rondror/akma327/GPCRContacts/data/fingerprint/general/gpcr_tree.json"
homology_model_path = "/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/homology-models"
classA_gpcr_path = "/scratch/PI/rondror/akma327/GPCRContacts/data/atomic-static-contacts/classA-gpcr-pdbs"

def get_uniprot_clusters():
	"""
		Return clusters of uniprot id that are cluster at the leaf node layer
	"""
	json_str = ""
	for line in open(gpcr_tree_path, 'r'):
		json_str += line 
	gpcr_tree = json.loads(json_str)

	uniprot_clusters = []
	
	for a in gpcr_tree["children"]:
		for b in a["children"]:
			cluster = []
			for uniprot in b["children"]:
				uniprot_id = (str(uniprot['name']) + "_human").upper()
				cluster.append(uniprot_id)
			uniprot_clusters.append(cluster)

	return uniprot_clusters


def driver(OUTPUT_DIR):
	uniprot_clusters = get_uniprot_clusters()
	for i, cluster in enumerate(uniprot_clusters):
		print("cluster_" + str(i))
		fp = []
		for uniprot in cluster:
			uni_id = uniprot.split("_")[0]
			hom_model_paths = glob.glob(homology_model_path + "/" + uni_id + "*/*table.txt")
			classA_paths = glob.glob(classA_gpcr_path + "/" + uni_id + "*/*table.txt")
			all_paths = hom_model_paths + classA_paths 
			fp += all_paths

		### Run gpcr_contact_fingerprint.py
		output_path = OUTPUT_DIR + "/" + "cluster_" + str(i)
		cmd_str = "python gpcr_contact_fingerprint.py " + output_path
		for p in fp:
			cmd_str += " " + p
		os.chdir("/scratch/PI/rondror/akma327/GPCRContacts/src/fingerprint")
		os.system(cmd_str)




if __name__ == "__main__":
	if(len(sys.argv) < K_MIN_ARG):
		print(USAGE_STR)
		exit(0)

	(OUTPUT_DIR) = (sys.argv[1])
	driver(OUTPUT_DIR)



