# Author: Anthony Kai Kwang Ma
# Date: 06/06/17
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# match_pdb_to_ligand.py

import sys
import glob

USAGE_STR = """
# Purpose
# Find mapping between pdb to ligand 


# Usage 
# python match_pdb_to_ligand.py 

"""

CLASSA_GPCR_PATH="/scratch/PI/rondror/akma327/GPCRContacts/data/structures/classA-gpcr-pdbs"
PDB_TO_LIGAND_PATH = "/scratch/PI/rondror/akma327/GPCRContacts/data/structures/classA-gpcr-pdbs/pdb_to_ligand.tsv"


def pdb_to_ligand():
	"""
		Create mapping from pdb to ligand 
	"""
	pdb_to_ligand_dict = {}
	f = open(PDB_TO_LIGAND_PATH, 'r')
	for line in f:
		uniprot, class_code, pdb, ligand = line.strip().split("\t")
		pdb_to_ligand_dict[pdb] = ligand

	return pdb_to_ligand_dict


def find_missing_ligands():
	"""
		Find ligands corresponding to pdb that are not yet assigned one
	"""


	pdb_to_ligand_dict = pdb_to_ligand()

	gpcr_list = glob.glob(CLASSA_GPCR_PATH + "/*pdb")
	for gpcr_path in gpcr_list:
		uniprot_pdb = gpcr_path.split("/")[-1].strip(".pdb").split("_")
		uniprot = (uniprot_pdb[0] + "_" + uniprot_pdb[1]).lower()
		pdb = uniprot_pdb[2]
		if (pdb not in pdb_to_ligand_dict):
			print(uniprot + "\t" + pdb)




if __name__ == "__main__":
	find_missing_ligands()

