# Author: Anthony Kai Kwang Ma
# Date: 06/06/17
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# prep_classA_gpcrs.py


import sys
import json
# import pymol
# from pymol import cmd, stored

USAGE_STR = """
# Purpose
# Scrape the full data set of class A gpcrs with specified chain from gpcrdb


# Usage 
# python prep_classA_gpcrs.py 

"""

gpcr_structures_file = "/Users/anthony/Desktop/dror/dynamic-networks1/GPCRContacts/src/prep-structures/gpcr_structures.json"
output_path="/Users/anthony/Desktop/dror/dynamic-networks1/GPCRContacts/data/structures_v2/classA-gpcr-pdbs"

def split_chains(selection='(all)', prefix=None):
    '''
DESCRIPTION

    Create a single object for each chain in selection

SEE ALSO

    split_states, http://pymolwiki.org/index.php/Split_object
    '''
    count = 0
    models = cmd.get_object_list('(' + selection + ')')
    for model in models:
        for chain in cmd.get_chains('(%s) and model %s' % (selection, model)):
            if chain == '':
                chain = "''"
            count += 1
            if not prefix:
                name = '%s_%s' % (model, chain)
            else:
                name = '%s%04d' % (prefix, count)
            cmd.create(name, '(%s) and model %s and chain %s' % (selection, model, chain))
        cmd.disable(model)

cmd.extend('split_chains', split_chains)


def get_gpcr_list():
	"""
		Acquire list of (uniprot, pdb, chain) tuples
	"""

	gpcr_list = []
	json_str = ""
	for line in open(gpcr_structures_file, 'r'):
		json_str += line

	gpcr_structures = json.loads(json_str)

	for entry in gpcr_structures:
		uniprot, pdb, chain = str(entry["protein"]).upper(), str(entry["pdb_code"]).upper(), str(entry["preferred_chain"]).split(",")[0]
		gpcr_list.append((uniprot, pdb, chain))

	return gpcr_list


def process_pdbs():
	gpcr_list = get_gpcr_list()

	for i, gpcr_tup in enumerate(gpcr_list):
		# if(i > 10): break
		uniprot, pdb, chain = gpcr_tup

		cmd.fetch(pdb)
		cmd.h_add("all")
		split_chains()
		chains = cmd.get_object_list()[1:]

		for pdb_c in chains:
			if(chain == pdb_c.split("_")[1]):
				cmd.save(output_path + "/" + uniprot + "_" + pdb_c + ".pdb", pdb_c)

		cmd.reinitialize()


process_pdbs()

