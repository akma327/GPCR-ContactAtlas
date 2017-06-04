# utils.py

PDB_TO_UNIPROT_TABLE_PATH = "/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/GPCR_PDB_List.txt"
GPCRDB_TABLE_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/All_species_gpcrdb_numbers_strOnly.txt"


def flipGpcrdbs(gpcrdb1, gpcrdb2):
       if("LIG" in gpcrdb1): return False
       elif("LIG" in gpcrdb2): return True

       tm1, index1 = map(int, gpcrdb1.split("x"))
       tm2, index2 = map(int, gpcrdb2.split("x"))

       if(tm1 < tm2): 
              return False
       elif(tm1 == tm2):
              if(index1 < index2):
                     return False
              else:
                     return True
       else:
              return True

def orderGpcrdbsAndResiAndAtom(gpcrdb1, gpcrdb2, resi1, resi2, atom1, atom2):
       """
              Return gpcrdb1 < gpcrdb2 ordered
       """

       if("LIG" in gpcrdb1): return (gpcrdb1, gpcrdb2, resi1, resi2, atom1, atom2)
       elif("LIG" in gpcrdb2): return (gpcrdb2, gpcrdb1, resi2, resi1, atom2, atom1)

       tm1, index1 = map(int, gpcrdb1.split("x"))
       tm2, index2 = map(int, gpcrdb2.split("x"))
       if(tm1 < tm2): 
              return (gpcrdb1, gpcrdb2, resi1, resi2, atom1, atom2)
       elif(tm1 == tm2): 
              if(index1 < index2): 
                     return (gpcrdb1, gpcrdb2, resi1, resi2, atom1, atom2)
              else: 
                     return (gpcrdb2, gpcrdb1, resi2, resi1, atom2, atom1)
       else:
              return (gpcrdb2, gpcrdb1, resi2, resi1, atom2, atom1)

def orderGpcrdbsAndResi(gpcrdb1, gpcrdb2, resi1, resi2):
       """
              Return gpcrdb1 < gpcrdb2 ordered
       """

       if("LIG" in gpcrdb1): return (gpcrdb1, gpcrdb2, resi1, resi2)
       elif("LIG" in gpcrdb2): return (gpcrdb2, gpcrdb1, resi2, resi1)

       tm1, index1 = map(int, gpcrdb1.split("x"))
       tm2, index2 = map(int, gpcrdb2.split("x"))
       if(tm1 < tm2): 
              return (gpcrdb1, gpcrdb2, resi1, resi2)
       elif(tm1 == tm2): 
              if(index1 < index2): 
                     return (gpcrdb1, gpcrdb2, resi1, resi2)
              else: 
                     return (gpcrdb2, gpcrdb1, resi2, resi1)
       else:
              return (gpcrdb2, gpcrdb1, resi2, resi1)

def orderGpcrdbs(gpcrdb1, gpcrdb2):
	"""
		Return gpcrdb1 < gpcrdb2 ordered
	"""

	if("LIG" in gpcrdb1): return (gpcrdb1, gpcrdb2)
	elif("LIG" in gpcrdb2): return (gpcrdb2, gpcrdb1)

	tm1, index1 = map(int, gpcrdb1.split("x"))
	tm2, index2 = map(int, gpcrdb2.split("x"))
	if(tm1 < tm2): 
		return (gpcrdb1, gpcrdb2)
	elif(tm1 == tm2): 
		if(index1 < index2): 
			return (gpcrdb1, gpcrdb2)
		else: 
			return (gpcrdb2, gpcrdb1)
	else:
		return (gpcrdb2, gpcrdb1)


def ligOrInTM(gpcrdb):
	if("LIGxLIG" == gpcrdb): return True
	tm, index = map(int, gpcrdb.split("x"))
	if(tm in range(1,9)): return True
	return False



# Rename amino acids to common name
def fixAminoAcidNames(key):
	key = key.replace("HSD", "HIS")
	key = key.replace("HSE", "HIS")
	key = key.replace("HSP", "HIS")
	key = key.replace("HIE", "HIS")
	key = key.replace("HIP", "HIS")
	key = key.replace("HID", "HIS")
	key = key.replace("GLH", "GLU")
	key = key.replace("ASH", "ASP")
	key = key.replace("CYP", "CYS")
	key = key.replace("CYX", "CYS")
	return key


# Retrive gpcrdb from dictionary for specified residue. Return "-" if not found
def getGPCRDB(res, GPCRDB_DICT):
	if("LIG" in res):
		return "LIGxLIG"
	res = fixAminoAcidNames(res)
	if(res not in GPCRDB_DICT):
		print(res + " not found.")
		return "None"
	return GPCRDB_DICT[res]


# Retrieve Uniprot Code for the PDB_CODE from pdb_to_uniprot_table_path
def getUniprotCode(PDB_CODE):
	f = open(PDB_TO_UNIPROT_TABLE_PATH, 'r')
	for line in f:
		if(line == "\n"): continue 
		l_info = line.split("\t")
		uniprot_code, pdb = l_info[0].strip(), l_info[2].strip()
		if(PDB_CODE.upper() == pdb.upper()): return uniprot_code.upper()
	print("PDB_CODE Not Found in PDB To Uniprot Table")
	exit(1)


# Given uniprot code reads through GPCRDB_TABLE_PATH to generate the amino acid
# to gpcrdb number table. 
# Output {"ASP112": "1x50", "ARG116":"2x45"}
def genGpcrdbDict(UNIPROT_CODE):
	GPCRDB_DICT = {}
	f = open(GPCRDB_TABLE_PATH, 'r')
	for line in f: 
		l_info = line.split("\t")
		uniprot, resnum, resname, gpcrdb = l_info[0].strip(), l_info[1].strip(), l_info[2].strip(), l_info[4].strip()
		if(uniprot.upper() == UNIPROT_CODE.upper()):
			key = resname.upper() + resnum 
			GPCRDB_DICT[key] = gpcrdb
	return GPCRDB_DICT


# Generates the residue to gpcrdb table for given pdb
def genResidueToGpcrdbTable(PDB_CODE):
	UNIPROT_CODE = getUniprotCode(PDB_CODE)
	GPCRDB_DICT = genGpcrdbDict(UNIPROT_CODE)
	return GPCRDB_DICT


node_dict =  {"tracks":[
    {"trackName":"Secondary Structure",
     "trackProperties":[
       {"nodeName":"1x28", "color":"#1500D6", "size":1},
       {"nodeName":"1x29", "color":"#1500D6", "size":1},
       {"nodeName":"1x30", "color":"#1500D6", "size":1},
       {"nodeName":"1x31", "color":"#1500D6", "size":1},
       {"nodeName":"1x32", "color":"#1500D6", "size":1},
       {"nodeName":"1x33", "color":"#1500D6", "size":1},
       {"nodeName":"1x34", "color":"#1500D6", "size":1},
       {"nodeName":"1x35", "color":"#1500D6", "size":1},
       {"nodeName":"1x36", "color":"#1500D6", "size":1},
       {"nodeName":"1x37", "color":"#1500D6", "size":1},
       {"nodeName":"1x38", "color":"#1500D6", "size":1},
       {"nodeName":"1x39", "color":"#1500D6", "size":1},
       {"nodeName":"1x40", "color":"#1500D6", "size":1},
       {"nodeName":"1x42", "color":"#1500D6", "size":1},
       {"nodeName":"1x43", "color":"#1500D6", "size":1},
       {"nodeName":"1x44", "color":"#1500D6", "size":1},
       {"nodeName":"1x45", "color":"#1500D6", "size":1},
       {"nodeName":"1x46", "color":"#1500D6", "size":1},
       {"nodeName":"1x47", "color":"#1500D6", "size":1},
       {"nodeName":"1x48", "color":"#1500D6", "size":1},
       {"nodeName":"1x49", "color":"#1500D6", "size":1},
       {"nodeName":"1x50", "color":"#1500D6", "size":1},
       {"nodeName":"1x51", "color":"#1500D6", "size":1},
       {"nodeName":"1x52", "color":"#1500D6", "size":1},
       {"nodeName":"1x53", "color":"#1500D6", "size":1},
       {"nodeName":"1x54", "color":"#1500D6", "size":1},
       {"nodeName":"1x55", "color":"#1500D6", "size":1},
       {"nodeName":"1x56", "color":"#1500D6", "size":1},
       {"nodeName":"1x57", "color":"#1500D6", "size":1},
       {"nodeName":"1x58", "color":"#1500D6", "size":1},
       {"nodeName":"1x59", "color":"#1500D6", "size":1},
       {"nodeName":"1x60", "color":"#1500D6", "size":1},
       {"nodeName":"2x38", "color":"#003D97", "size":1},
       {"nodeName":"2x39", "color":"#003D97", "size":1},
       {"nodeName":"2x40", "color":"#003D97", "size":1},
       {"nodeName":"2x41", "color":"#003D97", "size":1},
       {"nodeName":"2x42", "color":"#003D97", "size":1},
       {"nodeName":"2x43", "color":"#003D97", "size":1},
       {"nodeName":"2x44", "color":"#003D97", "size":1},
       {"nodeName":"2x45", "color":"#003D97", "size":1},
       {"nodeName":"2x46", "color":"#003D97", "size":1},
       {"nodeName":"2x47", "color":"#003D97", "size":1},
       {"nodeName":"2x48", "color":"#003D97", "size":1},
       {"nodeName":"2x49", "color":"#003D97", "size":1},
       {"nodeName":"2x50", "color":"#003D97", "size":1},
       {"nodeName":"2x51", "color":"#003D97", "size":1},
       {"nodeName":"2x52", "color":"#003D97", "size":1},
       {"nodeName":"2x53", "color":"#003D97", "size":1},
       {"nodeName":"2x54", "color":"#003D97", "size":1},
       {"nodeName":"2x55", "color":"#003D97", "size":1},
       {"nodeName":"2x56", "color":"#003D97", "size":1},
       {"nodeName":"2x57", "color":"#003D97", "size":1},
       {"nodeName":"2x58", "color":"#003D97", "size":1},
       {"nodeName":"2x59", "color":"#003D97", "size":1},
       {"nodeName":"2x60", "color":"#003D97", "size":1},
       {"nodeName":"2x61", "color":"#003D97", "size":1},
       {"nodeName":"2x62", "color":"#003D97", "size":1},
       {"nodeName":"2x63", "color":"#003D97", "size":1},
       {"nodeName":"2x64", "color":"#003D97", "size":1},
       {"nodeName":"2x65", "color":"#003D97", "size":1},
       {"nodeName":"2x66", "color":"#003D97", "size":1},
       {"nodeName":"3x22", "color":"#00E600", "size":1},
       {"nodeName":"3x23", "color":"#00E600", "size":1},
       {"nodeName":"3x24", "color":"#00E600", "size":1},
       {"nodeName":"3x25", "color":"#00E600", "size":1},
       {"nodeName":"3x26", "color":"#00E600", "size":1},
       {"nodeName":"3x27", "color":"#00E600", "size":1},
       {"nodeName":"3x28", "color":"#00E600", "size":1},
       {"nodeName":"3x29", "color":"#00E600", "size":1},
       {"nodeName":"3x30", "color":"#00E600", "size":1},
       {"nodeName":"3x31", "color":"#00E600", "size":1},
       {"nodeName":"3x32", "color":"#00E600", "size":1},
       {"nodeName":"3x33", "color":"#00E600", "size":1},
       {"nodeName":"3x34", "color":"#00E600", "size":1},
       {"nodeName":"3x35", "color":"#00E600", "size":1},
       {"nodeName":"3x36", "color":"#00E600", "size":1},
       {"nodeName":"3x37", "color":"#00E600", "size":1},
       {"nodeName":"3x38", "color":"#00E600", "size":1},
       {"nodeName":"3x39", "color":"#00E600", "size":1},
       {"nodeName":"3x40", "color":"#00E600", "size":1},
       {"nodeName":"3x41", "color":"#00E600", "size":1},
       {"nodeName":"3x42", "color":"#00E600", "size":1},
       {"nodeName":"3x43", "color":"#00E600", "size":1},
       {"nodeName":"3x44", "color":"#00E600", "size":1},
       {"nodeName":"3x45", "color":"#00E600", "size":1},
       {"nodeName":"3x46", "color":"#00E600", "size":1},
       {"nodeName":"3x47", "color":"#00E600", "size":1},
       {"nodeName":"3x48", "color":"#00E600", "size":1},
       {"nodeName":"3x49", "color":"#00E600", "size":1},
       {"nodeName":"3x50", "color":"#00E600", "size":1},
       {"nodeName":"3x51", "color":"#00E600", "size":1},
       {"nodeName":"3x52", "color":"#00E600", "size":1},
       {"nodeName":"3x53", "color":"#00E600", "size":1},
       {"nodeName":"3x54", "color":"#00E600", "size":1},
       {"nodeName":"3x55", "color":"#00E600", "size":1},
       {"nodeName":"3x56", "color":"#00E600", "size":1},
       {"nodeName":"4x41", "color":"#00E600", "size":1},
       {"nodeName":"4x42", "color":"#00E600", "size":1},
       {"nodeName":"4x43", "color":"#00E600", "size":1},
       {"nodeName":"4x44", "color":"#00E600", "size":1},
       {"nodeName":"4x45", "color":"#00E600", "size":1},
       {"nodeName":"4x46", "color":"#00E600", "size":1},
       {"nodeName":"4x47", "color":"#00E600", "size":1},
       {"nodeName":"4x48", "color":"#00E600", "size":1},
       {"nodeName":"4x49", "color":"#00E600", "size":1},
       {"nodeName":"4x50", "color":"#00E600", "size":1},
       {"nodeName":"4x51", "color":"#00E600", "size":1},
       {"nodeName":"4x52", "color":"#00E600", "size":1},
       {"nodeName":"4x53", "color":"#00E600", "size":1},
       {"nodeName":"4x54", "color":"#00E600", "size":1},
       {"nodeName":"4x55", "color":"#00E600", "size":1},
       {"nodeName":"4x56", "color":"#00E600", "size":1},
       {"nodeName":"4x57", "color":"#00E600", "size":1},
       {"nodeName":"4x58", "color":"#00E600", "size":1},
       {"nodeName":"4x59", "color":"#00E600", "size":1},
       {"nodeName":"4x60", "color":"#00E600", "size":1},
       {"nodeName":"4x61", "color":"#00E600", "size":1},
       {"nodeName":"4x62", "color":"#00E600", "size":1},
       {"nodeName":"4x63", "color":"#00E600", "size":1},
       {"nodeName":"5x32", "color":"#FEE200", "size":1},
       {"nodeName":"5x33", "color":"#FEE200", "size":1},
       {"nodeName":"5x34", "color":"#FEE200", "size":1},
       {"nodeName":"5x35", "color":"#FEE200", "size":1},
       {"nodeName":"5x36", "color":"#FEE200", "size":1},
       {"nodeName":"5x37", "color":"#FEE200", "size":1},
       {"nodeName":"5x38", "color":"#FEE200", "size":1},
       {"nodeName":"5x39", "color":"#FEE200", "size":1},
       {"nodeName":"5x40", "color":"#FEE200", "size":1},
       {"nodeName":"5x41", "color":"#FEE200", "size":1},
       {"nodeName":"5x42", "color":"#FEE200", "size":1},
       {"nodeName":"5x43", "color":"#FEE200", "size":1},
       {"nodeName":"5x44", "color":"#FEE200", "size":1},
       {"nodeName":"5x45", "color":"#FEE200", "size":1},
       {"nodeName":"5x46", "color":"#FEE200", "size":1},
       {"nodeName":"5x47", "color":"#FEE200", "size":1},
       {"nodeName":"5x48", "color":"#FEE200", "size":1},
       {"nodeName":"5x49", "color":"#FEE200", "size":1},
       {"nodeName":"5x50", "color":"#FEE200", "size":1},
       {"nodeName":"5x51", "color":"#FEE200", "size":1},
       {"nodeName":"5x52", "color":"#FEE200", "size":1},
       {"nodeName":"5x53", "color":"#FEE200", "size":1},
       {"nodeName":"5x54", "color":"#FEE200", "size":1},
       {"nodeName":"5x55", "color":"#FEE200", "size":1},
       {"nodeName":"5x56", "color":"#FEE200", "size":1},
       {"nodeName":"5x57", "color":"#FEE200", "size":1},
       {"nodeName":"5x58", "color":"#FEE200", "size":1},
       {"nodeName":"5x59", "color":"#FEE200", "size":1},
       {"nodeName":"5x60", "color":"#FEE200", "size":1},
       {"nodeName":"5x61", "color":"#FEE200", "size":1},
       {"nodeName":"5x62", "color":"#FEE200", "size":1},
       {"nodeName":"5x63", "color":"#FEE200", "size":1},
       {"nodeName":"5x64", "color":"#FEE200", "size":1},
       {"nodeName":"5x461","color":"#FEE200", "size":1},
       {"nodeName":"6x24", "color":"#FF9000", "size":1},
       {"nodeName":"6x25", "color":"#FF9000", "size":1},
       {"nodeName":"6x26", "color":"#FF9000", "size":1},
       {"nodeName":"6x27", "color":"#FF9000", "size":1},
       {"nodeName":"6x28", "color":"#FF9000", "size":1},
       {"nodeName":"6x29", "color":"#FF9000", "size":1},
       {"nodeName":"6x30", "color":"#FF9000", "size":1},
       {"nodeName":"6x31", "color":"#FF9000", "size":1},
       {"nodeName":"6x32", "color":"#FF9000", "size":1},
       {"nodeName":"6x33", "color":"#FF9000", "size":1},
       {"nodeName":"6x34", "color":"#FF9000", "size":1},
       {"nodeName":"6x35", "color":"#FF9000", "size":1},
       {"nodeName":"6x36", "color":"#FF9000", "size":1},
       {"nodeName":"6x37", "color":"#FF9000", "size":1},
       {"nodeName":"6x38", "color":"#FF9000", "size":1},
       {"nodeName":"6x39", "color":"#FF9000", "size":1},
       {"nodeName":"6x40", "color":"#FF9000", "size":1},
       {"nodeName":"6x41", "color":"#FF9000", "size":1},
       {"nodeName":"6x42", "color":"#FF9000", "size":1},
       {"nodeName":"6x43", "color":"#FF9000", "size":1},
       {"nodeName":"6x44", "color":"#FF9000", "size":1},
       {"nodeName":"6x45", "color":"#FF9000", "size":1},
       {"nodeName":"6x46", "color":"#FF9000", "size":1},
       {"nodeName":"6x47", "color":"#FF9000", "size":1},
       {"nodeName":"6x48", "color":"#FF9000", "size":1},
       {"nodeName":"6x49", "color":"#FF9000", "size":1},
       {"nodeName":"6x50", "color":"#FF9000", "size":1},
       {"nodeName":"6x51", "color":"#FF9000", "size":1},
       {"nodeName":"6x52", "color":"#FF9000", "size":1},
       {"nodeName":"6x53", "color":"#FF9000", "size":1},
       {"nodeName":"6x54", "color":"#FF9000", "size":1},
       {"nodeName":"6x55", "color":"#FF9000", "size":1},
       {"nodeName":"6x56", "color":"#FF9000", "size":1},
       {"nodeName":"6x57", "color":"#FF9000", "size":1},
       {"nodeName":"6x58", "color":"#FF9000", "size":1},
       {"nodeName":"6x59", "color":"#FF9000", "size":1},
       {"nodeName":"6x60", "color":"#FF9000", "size":1},
       {"nodeName":"6x61", "color":"#FF9000", "size":1},
       {"nodeName":"7x28", "color":"#FF3B00", "size":1},
       {"nodeName":"7x29", "color":"#FF3B00", "size":1},
       {"nodeName":"7x30", "color":"#FF3B00", "size":1},
       {"nodeName":"7x31", "color":"#FF3B00", "size":1},
       {"nodeName":"7x32", "color":"#FF3B00", "size":1},
       {"nodeName":"7x33", "color":"#FF3B00", "size":1},
       {"nodeName":"7x34", "color":"#FF3B00", "size":1},
       {"nodeName":"7x35", "color":"#FF3B00", "size":1},
       {"nodeName":"7x36", "color":"#FF3B00", "size":1},
       {"nodeName":"7x37", "color":"#FF3B00", "size":1},
       {"nodeName":"7x38", "color":"#FF3B00", "size":1},
       {"nodeName":"7x39", "color":"#FF3B00", "size":1},
       {"nodeName":"7x40", "color":"#FF3B00", "size":1},
       {"nodeName":"7x41", "color":"#FF3B00", "size":1},
       {"nodeName":"7x42", "color":"#FF3B00", "size":1},
       {"nodeName":"7x43", "color":"#FF3B00", "size":1},
       {"nodeName":"7x44", "color":"#FF3B00", "size":1},
       {"nodeName":"7x45", "color":"#FF3B00", "size":1},
       {"nodeName":"7x46", "color":"#FF3B00", "size":1},
       {"nodeName":"7x47", "color":"#FF3B00", "size":1},
       {"nodeName":"7x48", "color":"#FF3B00", "size":1},
       {"nodeName":"7x49", "color":"#FF3B00", "size":1},
       {"nodeName":"7x50", "color":"#FF3B00", "size":1},
       {"nodeName":"7x51", "color":"#FF3B00", "size":1},
       {"nodeName":"7x52", "color":"#FF3B00", "size":1},
       {"nodeName":"7x53", "color":"#FF3B00", "size":1},
       {"nodeName":"7x54", "color":"#FF3B00", "size":1},
       {"nodeName":"7x55", "color":"#FF3B00", "size":1},
       {"nodeName":"7x56", "color":"#FF3B00", "size":1},
       {"nodeName":"7x57", "color":"#FF3B00", "size":1},
       {"nodeName":"8x48", "color":"#FF0000", "size":1},
       {"nodeName":"8x49", "color":"#FF0000", "size":1},
       {"nodeName":"8x50", "color":"#FF0000", "size":1},
       {"nodeName":"8x51", "color":"#FF0000", "size":1},
       {"nodeName":"8x52", "color":"#FF0000", "size":1},
       {"nodeName":"8x53", "color":"#FF0000", "size":1},
       {"nodeName":"8x54", "color":"#FF0000", "size":1}
      ]
    }
  ],

  "trees":[
    {"treeName":"Secondary Structure",
     "treePaths":[
       "1.1x28",
       "1.1x29",
       "1.1x30",
       "1.1x31",
       "1.1x32",
       "1.1x33",
       "1.1x34",
       "1.1x35",
       "1.1x36",
       "1.1x37",
       "1.1x38",
       "1.1x39",
       "1.1x40",
       "1.1x42",
       "1.1x43",
       "1.1x44",
       "1.1x45",
       "1.1x46",
       "1.1x47",
       "1.1x48",
       "1.1x49",
       "1.1x50",
       "1.1x51",
       "1.1x52",
       "1.1x53",
       "1.1x54",
       "1.1x55",
       "1.1x56",
       "1.1x57",
       "1.1x58",
       "1.1x59",
       "1.1x60",
       "2.2x38",
       "2.2x39",
       "2.2x40",
       "2.2x41",
       "2.2x42",
       "2.2x43",
       "2.2x44",
       "2.2x45",
       "2.2x46",
       "2.2x47",
       "2.2x48",
       "2.2x49",
       "2.2x50",
       "2.2x51",
       "2.2x52",
       "2.2x53",
       "2.2x54",
       "2.2x55",
       "2.2x56",
       "2.2x57",
       "2.2x58",
       "2.2x59",
       "2.2x60",
       "2.2x61",
       "2.2x62",
       "2.2x63",
       "2.2x64",
       "2.2x65",
       "2.2x66",
       "3.3x22",
       "3.3x23",
       "3.3x24",
       "3.3x25",
       "3.3x26",
       "3.3x27",
       "3.3x28",
       "3.3x29",
       "3.3x30",
       "3.3x31",
       "3.3x32",
       "3.3x33",
       "3.3x34",
       "3.3x35",
       "3.3x36",
       "3.3x37",
       "3.3x38",
       "3.3x39",
       "3.3x40",
       "3.3x41",
       "3.3x42",
       "3.3x43",
       "3.3x44",
       "3.3x45",
       "3.3x46",
       "3.3x47",
       "3.3x48",
       "3.3x49",
       "3.3x50",
       "3.3x51",
       "3.3x52",
       "3.3x53",
       "3.3x54",
       "3.3x55",
       "3.3x56",
       "4.4x41",
       "4.4x42",
       "4.4x43",
       "4.4x44",
       "4.4x45",
       "4.4x46",
       "4.4x47",
       "4.4x48",
       "4.4x49",
       "4.4x50",
       "4.4x51",
       "4.4x52",
       "4.4x53",
       "4.4x54",
       "4.4x55",
       "4.4x56",
       "4.4x57",
       "4.4x58",
       "4.4x59",
       "4.4x60",
       "4.4x61",
       "4.4x62",
       "4.4x63",
       "5.5x32",
       "5.5x33",
       "5.5x34",
       "5.5x35",
       "5.5x36",
       "5.5x37",
       "5.5x38",
       "5.5x39",
       "5.5x40",
       "5.5x41",
       "5.5x42",
       "5.5x43",
       "5.5x44",
       "5.5x45",
       "5.5x46",
       "5.5x47",
       "5.5x48",
       "5.5x49",
       "5.5x50",
       "5.5x51",
       "5.5x52",
       "5.5x53",
       "5.5x54",
       "5.5x55",
       "5.5x56",
       "5.5x57",
       "5.5x58",
       "5.5x59",
       "5.5x60",
       "5.5x61",
       "5.5x62",
       "5.5x63",
       "5.5x64",
       "5.5x461",
       "6.6x24",
       "6.6x25",
       "6.6x26",
       "6.6x27",
       "6.6x28",
       "6.6x29",
       "6.6x30",
       "6.6x31",
       "6.6x32",
       "6.6x33",
       "6.6x34",
       "6.6x35",
       "6.6x36",
       "6.6x37",
       "6.6x38",
       "6.6x39",
       "6.6x40",
       "6.6x41",
       "6.6x42",
       "6.6x43",
       "6.6x44",
       "6.6x45",
       "6.6x46",
       "6.6x47",
       "6.6x48",
       "6.6x49",
       "6.6x50",
       "6.6x51",
       "6.6x52",
       "6.6x53",
       "6.6x54",
       "6.6x55",
       "6.6x56",
       "6.6x57",
       "6.6x58",
       "6.6x59",
       "6.6x60",
       "6.6x61",
       "7.7x28",
       "7.7x29",
       "7.7x30",
       "7.7x31",
       "7.7x32",
       "7.7x33",
       "7.7x34",
       "7.7x35",
       "7.7x36",
       "7.7x37",
       "7.7x38",
       "7.7x39",
       "7.7x40",
       "7.7x41",
       "7.7x42",
       "7.7x43",
       "7.7x44",
       "7.7x45",
       "7.7x46",
       "7.7x47",
       "7.7x48",
       "7.7x49",
       "7.7x50",
       "7.7x51",
       "7.7x52",
       "7.7x53",
       "7.7x54",
       "7.7x55",
       "7.7x56",
       "7.7x57",
       "8.8x48",
       "8.8x49",
       "8.8x50",
       "8.8x51",
       "8.8x52",
       "8.8x53",
       "8.8x54"
      ]
    }
  ],
  "defaults":{
    "edgeColor":"rgba(100,100,100,100)",
    "edgeWidth":1
  }
}
