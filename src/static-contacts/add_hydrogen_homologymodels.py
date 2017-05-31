# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# Date: 04/29/17
# add_hydrogen_homologymodels.py

import pymol
from pymol import cmd, stored


import glob 


USAGE_STR = """
# Purpose
# Generate .sh script to call StaticInteractionCalculator.py upon every available
# homology model or existing crystal structure. 

# Usage 
# python add_hydrogen_homologymodels.py 

"""


HOMOLOGY_MODEL_PATH="/Users/anthony/Desktop/dror/dynamic-networks1/GPCRContacts/data/homology_models"
HADDED_HOMOLOGY_MODEL_PATH="/Users/anthony/Desktop/dror/dynamic-networks1/GPCRContacts/data/hadded-homology-models"

def add_hydrogens():
	model_paths = glob.glob(HOMOLOGY_MODEL_PATH + "/*/*pdb")
	for i, mp in enumerate(model_paths):
		# if(i > 3): break
		label = "_".join(mp.split("/")[-1].split("_")[0:2]).upper()
		output_path = HADDED_HOMOLOGY_MODEL_PATH + "/" + label + ".pdb"

		### Retrieve PDB and add hydrogens
		cmd.load(mp, label)
		cmd.h_add("all")
		cmd.save(output_path, label)

		cmd.reinitialize()



add_hydrogens()