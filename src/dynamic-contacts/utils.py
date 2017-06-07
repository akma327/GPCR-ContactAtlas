# Author: Anthony Kai Kwang Ma
# Email: anthonyma27@gmail.com, akma327@stanford.edu
# utils.py 

import os 
import sys
import csv
import numpy as np

# Paths to important configuration tables
PDB_TO_UNIPROT_TABLE_PATH = "/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/GPCR_PDB_List.txt"
GPCRDB_TABLE_PATH="/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/All_species_gpcrdb_numbers_strOnly.txt"
GPCRDB_RESIDUE_FREQ_TABLE="/scratch/PI/rondror/akma327/DynamicNetworks/data/crystal-analysis/simulation-analysis/gpcrdb-freq-config/gpcrdb_residue_conservation.txt"

# Single letter to three letter representation of amino acids
aa1to3_dict = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}
aa3to1_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'ASN': 'N', 'GLN': 'Q', 'LYS': 'K', 'THR': 'T', 'PRO': 'P', 'HIS': 'H', 'PHE': 'F', 'ALA': 'A', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# Flag to interaction dictionary file names 
flag_to_dict_file = {'sb':'salt_bridge_result_dict.txt', 'pc':'pi_cation_result_dict.txt', 'ps':'pi_stacking_result_dict.txt', 'ts':'t_stacking_result_dict.txt',
'vdw':'vanderwaal_result_dict.txt', 'hbbb':'bb_hydrogen_bond_result_dict.txt', 'hbsb':'sb_hydrogen_bond_result_dict.txt', 'hbss':'ss_hydrogen_bond_result_dict.txt',
'rw':'residue_water_hydrogen_bond_result_dict.txt', 'wb':'water_bond_result_dict.txt', 'wb2':'extended_water_bond_result_dict.txt',
'hlb':'lb_hydrogen_bond_result_dict.txt', 'hls':'ls_hydrogen_bond_result_dict.txt', 'lw':'ligand_water_hydrogen_bond_result_dict.txt', 
'lwb':'ligand_water_bond_result_dict.txt', 'lwb2':'ligand_extended_water_bond_result_dict.txt'}

# Alphabetical list of 20 amino acids
ALPHABETICAL_AA_LIST = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

# Conservation Coloring Dictionary 
GRAY_SCALE_DICT={"gray32": "#525252", "gray33": "#545454", "gray30": "#4D4D4D", "gray31": "#4F4F4F", "gray36": "#5C5C5C", "gray37": "#5E5E5E", "gray34": "#575757", "gray35": "#595959", "gray38": "#616161", "gray39": "#636363", "gray8": "#141414", "gray9": "#171717", "gray2": "#050505", "gray3": "#080808", "gray1": "#030303", "gray6": "#0F0F0F", "gray7": "#121212", "gray4": "#0A0A0A", "gray5": "#0D0D0D", "gray21": "#363636", "gray20": "#333333", "gray23": "#3B3B3B", "gray22": "#383838", "gray25": "#404040", "gray24": "#3D3D3D", "gray27": "#454545", "gray26": "#424242", "gray29": "#4A4A4A", "gray28": "#474747", "gray18": "#2E2E2E", "gray19": "#303030", "gray14": "#242424", "gray15": "#262626", "gray16": "#292929", "gray17": "#2B2B2B", "gray10": "#1A1A1A", "gray11": "#1C1C1C", "gray12": "#1F1F1F", "gray13": "#212121", "gray94": "#F0F0F0", "gray95": "#F2F2F2", "gray96": "#F5F5F5", "gray97": "#F7F7F7", "gray90": "#E5E5E5", "gray91": "#E8E8E8", "gray92": "#EBEBEB", "gray93": "#EDEDED", "gray98": "#FAFAFA", "gray99": "#FCFCFC", "gray83": "#D4D4D4", "gray82": "#D1D1D1", "gray81": "#CFCFCF", "gray80": "#CCCCCC", "gray87": "#DEDEDE", "gray86": "#DBDBDB", "gray85": "#D9D9D9", "gray84": "#D6D6D6", "gray89": "#E3E3E3", "gray88": "#E0E0E0", "gray78": "#C7C7C7", "gray79": "#C9C9C9", "gray76": "#C2C2C2", "gray77": "#C4C4C4", "gray74": "#BDBDBD", "gray75": "#BFBFBF", "gray72": "#B8B8B8", "gray73": "#BABABA", "gray70": "#B3B3B3", "gray71": "#B5B5B5", "gray69": "#B0B0B0", "gray68": "#ADADAD", "gray65": "#A6A6A6", "gray64": "#A3A3A3", "gray67": "#ABABAB", "gray66": "#A8A8A8", "gray61": "#9C9C9C", "gray60": "#999999", "gray63": "#A1A1A1", "gray62": "#9E9E9E", "gray58": "#949494", "gray59": "#969696", "gray50": "#7F7F7F", "gray51": "#828282", "gray52": "#858585", "gray53": "#878787", "gray54": "#8A8A8A", "gray55": "#8C8C8C", "gray56": "#8F8F8F", "gray57": "#919191", "gray47": "#787878", "gray46": "#757575", "gray45": "#737373", "gray44": "#707070", "gray43": "#6E6E6E", "gray42": "#6B6B6B", "gray41": "#696969", "gray40": "#666666", "gray49": "#7D7D7D", "gray48": "#7A7A7A", "gray100": "#FFFFFF"}


# Bundle plot drawing color 
ORANGE="FFA500"
CYAN="00FFFF"
GREEN="#008000"


# COLOR CODING For different receptors
MOR = "00FFFF"
D4R = "FFFF00"
M2R = "FF00FF"
B2AR = "FFA500"
DOR = "#00FF00"

RECEPTOR_TO_COLOR = {"MOR": "00FFFF", "D4R": "FFFF00", "M2R" : "FF00FF", "B2AR" : "FFA500", "DOR" : "00FF00", "RHO": "800000", "US28": "#FF0000", "A2A": "#6600FF"}

LINE_COLOR_DICT = {"hbss":ORANGE , "hbbb": ORANGE, "hbsb": ORANGE, "wb": CYAN, "wb2": GREEN, "wb_combined" : B2AR}
LINE_COLOR_TO_TYPE_DICT = {ORANGE:"hb", CYAN:"wb", GREEN:"wb2"}


# Reads json file into json string 
def loadJsonString(POLAR_NETWORK_JSON1):
	f = open(POLAR_NETWORK_JSON1, 'r')
	json_str1 = ""
	for line in f:
		json_str1 += line
	return json_str1


# Find column index for EXP_REP_COLUMN_ID
def findColumnIndex(column_labels, EXP_REP_COLUMN_ID):
	for index, v in enumerate(column_labels):
		if (EXP_REP_COLUMN_ID in v):
			return index 
	print("EXP_REP_COLUMN_ID: " + EXP_REP_COLUMN_ID + " not found in INTERACTION_FOOTPRINT_TABLE.")
	exit(1)	


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
		return "------"
    return GPCRDB_DICT[res]

# Create directory if not exist
def createDirectory(OUTPUT_FILE):
	directory = os.path.dirname(OUTPUT_FILE)
	if not os.path.exists(directory):
		os.makedirs(directory)

# Generate write file descriptor 
def genWriteDescriptor(OUTPUT_FILE):
	createDirectory(OUTPUT_FILE)
	return open(OUTPUT_FILE, 'w')


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

def genGPCRDBToResidueTable(PDB_CODE):
    GPCRDB_DICT = genResidueToGpcrdbTable(PDB_CODE)
    gpcrdb_to_residue = {GPCRDB_DICT[resid] : resid for resid in GPCRDB_DICT}
    return gpcrdb_to_residue

# Input: Interaction Footprint Table and specified column index
# Function: Extract list of pairs of atoms and associated interaction type and frequency value 
# Output: List of tuples containing (interaction_row_label, frequency_value)
def processInteractionFootprintTable(INTERACTION_FOOTPRINT_TABLE, EXP_REP_COLUMN_ID):
	with open(INTERACTION_FOOTPRINT_TABLE, 'rU') as f: 
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	darr = np.array(d)
	column_labels = list(darr[0,:])
	column_index = findColumnIndex(column_labels, EXP_REP_COLUMN_ID)
	interaction_row_labels = list(darr[:,0])[1:]
	exp_rep_freq_vals = map(float, list(darr[:, column_index])[1:])
	row_label_freq = [(row_label, exp_rep_freq_vals[index]) for index, row_label in enumerate(interaction_row_labels)]
	return row_label_freq



# Parses "LYS384-NZ_1 -- GLU382-OE2_1" into ("LYS384", "GLU382")
def extractResPair(resAtomPair):
	resatom1, resatom2 = resAtomPair.split(" -- ")
	res1, res2 = resatom1.split("-")[0].strip(), resatom2.split("-")[0].strip()
	return res1, res2

# Parses "LYS384-NZ_1 -- GLU382-OE2_1" into (GPCRDB1, GPCRDB2)
def extractGpcrdbPair(resAtomPair, EXP_COND_PDB_DICT):
	res1, res2 = extractResPair(resAtomPair)
	return getGPCRDB(res1, EXP_COND_PDB_DICT), getGPCRDB(res2, EXP_COND_PDB_DICT)



NODE_COLOR = [
    {"name":"1.1x28", "color":"#1500D6"},
    {"name":"1.1x29", "color":"#1500D6"},
    {"name":"1.1x30", "color":"#1500D6"},
    {"name":"1.1x31", "color":"#1500D6"},
    {"name":"1.1x32", "color":"#1500D6"},
    {"name":"1.1x33", "color":"#1500D6"},
    {"name":"1.1x34", "color":"#1500D6"},
    {"name":"1.1x35", "color":"#1500D6"},
    {"name":"1.1x36", "color":"#1500D6"},
    {"name":"1.1x37", "color":"#1500D6"},
    {"name":"1.1x38", "color":"#1500D6"},
    {"name":"1.1x39", "color":"#1500D6"},
    {"name":"1.1x40", "color":"#1500D6"},
    {"name":"1.1x42", "color":"#1500D6"},
    {"name":"1.1x43", "color":"#1500D6"},
    {"name":"1.1x44", "color":"#1500D6"},
    {"name":"1.1x45", "color":"#1500D6"},
    {"name":"1.1x46", "color":"#1500D6"},
    {"name":"1.1x47", "color":"#1500D6"},
    {"name":"1.1x48", "color":"#1500D6"},
    {"name":"1.1x49", "color":"#1500D6"},
    {"name":"1.1x50", "color":"#1500D6"},
    {"name":"1.1x51", "color":"#1500D6"},
    {"name":"1.1x52", "color":"#1500D6"},
    {"name":"1.1x53", "color":"#1500D6"},
    {"name":"1.1x54", "color":"#1500D6"},
    {"name":"1.1x55", "color":"#1500D6"},
    {"name":"1.1x56", "color":"#1500D6"},
    {"name":"1.1x57", "color":"#1500D6"},
    {"name":"1.1x58", "color":"#1500D6"},
    {"name":"1.1x59", "color":"#1500D6"},
    {"name":"1.1x60", "color":"#1500D6"},
    {"name":"2.2x38", "color":"#003D97"},
    {"name":"2.2x39", "color":"#003D97"},
    {"name":"2.2x40", "color":"#003D97"},
    {"name":"2.2x41", "color":"#003D97"},
    {"name":"2.2x42", "color":"#003D97"},
    {"name":"2.2x43", "color":"#003D97"},
    {"name":"2.2x44", "color":"#003D97"},
    {"name":"2.2x45", "color":"#003D97"},
    {"name":"2.2x46", "color":"#003D97"},
    {"name":"2.2x47", "color":"#003D97"},
    {"name":"2.2x48", "color":"#003D97"},
    {"name":"2.2x49", "color":"#003D97"},
    {"name":"2.2x50", "color":"#003D97"},
    {"name":"2.2x51", "color":"#003D97"},
    {"name":"2.2x52", "color":"#003D97"},
    {"name":"2.2x53", "color":"#003D97"},
    {"name":"2.2x54", "color":"#003D97"},
    {"name":"2.2x55", "color":"#003D97"},
    {"name":"2.2x56", "color":"#003D97"},
    {"name":"2.2x57", "color":"#003D97"},
    {"name":"2.2x58", "color":"#003D97"},
    {"name":"2.2x59", "color":"#003D97"},
    {"name":"2.2x60", "color":"#003D97"},
    {"name":"2.2x61", "color":"#003D97"},
    {"name":"2.2x62", "color":"#003D97"},
    {"name":"2.2x63", "color":"#003D97"},
    {"name":"2.2x64", "color":"#003D97"},
    {"name":"2.2x65", "color":"#003D97"},
    {"name":"2.2x66", "color":"#003D97"},
    {"name":"3.3x22", "color":"#00E600"},
    {"name":"3.3x23", "color":"#00E600"},
    {"name":"3.3x24", "color":"#00E600"},
    {"name":"3.3x25", "color":"#00E600"},
    {"name":"3.3x26", "color":"#00E600"},
    {"name":"3.3x27", "color":"#00E600"},
    {"name":"3.3x28", "color":"#00E600"},
    {"name":"3.3x29", "color":"#00E600"},
    {"name":"3.3x30", "color":"#00E600"},
    {"name":"3.3x31", "color":"#00E600"},
    {"name":"3.3x32", "color":"#00E600"},
    {"name":"3.3x33", "color":"#00E600"},
    {"name":"3.3x34", "color":"#00E600"},
    {"name":"3.3x35", "color":"#00E600"},
    {"name":"3.3x36", "color":"#00E600"},
    {"name":"3.3x37", "color":"#00E600"},
    {"name":"3.3x38", "color":"#00E600"},
    {"name":"3.3x39", "color":"#00E600"},
    {"name":"3.3x40", "color":"#00E600"},
    {"name":"3.3x41", "color":"#00E600"},
    {"name":"3.3x42", "color":"#00E600"},
    {"name":"3.3x43", "color":"#00E600"},
    {"name":"3.3x44", "color":"#00E600"},
    {"name":"3.3x45", "color":"#00E600"},
    {"name":"3.3x46", "color":"#00E600"},
    {"name":"3.3x47", "color":"#00E600"},
    {"name":"3.3x48", "color":"#00E600"},
    {"name":"3.3x49", "color":"#00E600"},
    {"name":"3.3x50", "color":"#00E600"},
    {"name":"3.3x51", "color":"#00E600"},
    {"name":"3.3x52", "color":"#00E600"},
    {"name":"3.3x53", "color":"#00E600"},
    {"name":"3.3x54", "color":"#00E600"},
    {"name":"3.3x55", "color":"#00E600"},
    {"name":"3.3x56", "color":"#00E600"},
    {"name":"4.4x41", "color":"#00E600"},
    {"name":"4.4x42", "color":"#00E600"},
    {"name":"4.4x43", "color":"#00E600"},
    {"name":"4.4x44", "color":"#00E600"},
    {"name":"4.4x45", "color":"#00E600"},
    {"name":"4.4x46", "color":"#00E600"},
    {"name":"4.4x47", "color":"#00E600"},
    {"name":"4.4x48", "color":"#00E600"},
    {"name":"4.4x49", "color":"#00E600"},
    {"name":"4.4x50", "color":"#00E600"},
    {"name":"4.4x51", "color":"#00E600"},
    {"name":"4.4x52", "color":"#00E600"},
    {"name":"4.4x53", "color":"#00E600"},
    {"name":"4.4x54", "color":"#00E600"},
    {"name":"4.4x55", "color":"#00E600"},
    {"name":"4.4x56", "color":"#00E600"},
    {"name":"4.4x57", "color":"#00E600"},
    {"name":"4.4x58", "color":"#00E600"},
    {"name":"4.4x59", "color":"#00E600"},
    {"name":"4.4x60", "color":"#00E600"},
    {"name":"4.4x61", "color":"#00E600"},
    {"name":"4.4x62", "color":"#00E600"},
    {"name":"4.4x63", "color":"#00E600"},
    {"name":"5.5x32", "color":"#FEE200"},
    {"name":"5.5x33", "color":"#FEE200"},
    {"name":"5.5x34", "color":"#FEE200"},
    {"name":"5.5x35", "color":"#FEE200"},
    {"name":"5.5x36", "color":"#FEE200"},
    {"name":"5.5x37", "color":"#FEE200"},
    {"name":"5.5x38", "color":"#FEE200"},
    {"name":"5.5x39", "color":"#FEE200"},
    {"name":"5.5x40", "color":"#FEE200"},
    {"name":"5.5x41", "color":"#FEE200"},
    {"name":"5.5x42", "color":"#FEE200"},
    {"name":"5.5x43", "color":"#FEE200"},
    {"name":"5.5x44", "color":"#FEE200"},
    {"name":"5.5x45", "color":"#FEE200"},
    {"name":"5.5x46", "color":"#FEE200"},
    {"name":"5.5x47", "color":"#FEE200"},
    {"name":"5.5x48", "color":"#FEE200"},
    {"name":"5.5x49", "color":"#FEE200"},
    {"name":"5.5x50", "color":"#FEE200"},
    {"name":"5.5x51", "color":"#FEE200"},
    {"name":"5.5x52", "color":"#FEE200"},
    {"name":"5.5x53", "color":"#FEE200"},
    {"name":"5.5x54", "color":"#FEE200"},
    {"name":"5.5x55", "color":"#FEE200"},
    {"name":"5.5x56", "color":"#FEE200"},
    {"name":"5.5x57", "color":"#FEE200"},
    {"name":"5.5x58", "color":"#FEE200"},
    {"name":"5.5x59", "color":"#FEE200"},
    {"name":"5.5x60", "color":"#FEE200"},
    {"name":"5.5x61", "color":"#FEE200"},
    {"name":"5.5x62", "color":"#FEE200"},
    {"name":"5.5x63", "color":"#FEE200"},
    {"name":"5.5x64", "color":"#FEE200"},
    {"name":"5.5x461", "color":"#FEE200"},
    {"name":"6.6x24", "color":"#FF9000"},
    {"name":"6.6x25", "color":"#FF9000"},
    {"name":"6.6x26", "color":"#FF9000"},
    {"name":"6.6x27", "color":"#FF9000"},
    {"name":"6.6x28", "color":"#FF9000"},
    {"name":"6.6x29", "color":"#FF9000"},
    {"name":"6.6x30", "color":"#FF9000"},
    {"name":"6.6x31", "color":"#FF9000"},
    {"name":"6.6x32", "color":"#FF9000"},
    {"name":"6.6x33", "color":"#FF9000"},
    {"name":"6.6x34", "color":"#FF9000"},
    {"name":"6.6x35", "color":"#FF9000"},
    {"name":"6.6x36", "color":"#FF9000"},
    {"name":"6.6x37", "color":"#FF9000"},
    {"name":"6.6x38", "color":"#FF9000"},
    {"name":"6.6x39", "color":"#FF9000"},
    {"name":"6.6x40", "color":"#FF9000"},
    {"name":"6.6x41", "color":"#FF9000"},
    {"name":"6.6x42", "color":"#FF9000"},
    {"name":"6.6x43", "color":"#FF9000"},
    {"name":"6.6x44", "color":"#FF9000"},
    {"name":"6.6x45", "color":"#FF9000"},
    {"name":"6.6x46", "color":"#FF9000"},
    {"name":"6.6x47", "color":"#FF9000"},
    {"name":"6.6x48", "color":"#FF9000"},
    {"name":"6.6x49", "color":"#FF9000"},
    {"name":"6.6x50", "color":"#FF9000"},
    {"name":"6.6x51", "color":"#FF9000"},
    {"name":"6.6x52", "color":"#FF9000"},
    {"name":"6.6x53", "color":"#FF9000"},
    {"name":"6.6x54", "color":"#FF9000"},
    {"name":"6.6x55", "color":"#FF9000"},
    {"name":"6.6x56", "color":"#FF9000"},
    {"name":"6.6x57", "color":"#FF9000"},
    {"name":"6.6x58", "color":"#FF9000"},
    {"name":"6.6x59", "color":"#FF9000"},
    {"name":"6.6x60", "color":"#FF9000"},
    {"name":"6.6x61", "color":"#FF9000"},
    {"name":"7.7x28", "color":"#FF3B00"},
    {"name":"7.7x29", "color":"#FF3B00"},
    {"name":"7.7x30", "color":"#FF3B00"},
    {"name":"7.7x31", "color":"#FF3B00"},
    {"name":"7.7x32", "color":"#FF3B00"},
    {"name":"7.7x33", "color":"#FF3B00"},
    {"name":"7.7x34", "color":"#FF3B00"},
    {"name":"7.7x35", "color":"#FF3B00"},
    {"name":"7.7x36", "color":"#FF3B00"},
    {"name":"7.7x37", "color":"#FF3B00"},
    {"name":"7.7x38", "color":"#FF3B00"},
    {"name":"7.7x39", "color":"#FF3B00"},
    {"name":"7.7x40", "color":"#FF3B00"},
    {"name":"7.7x41", "color":"#FF3B00"},
    {"name":"7.7x42", "color":"#FF3B00"},
    {"name":"7.7x43", "color":"#FF3B00"},
    {"name":"7.7x44", "color":"#FF3B00"},
    {"name":"7.7x45", "color":"#FF3B00"},
    {"name":"7.7x46", "color":"#FF3B00"},
    {"name":"7.7x47", "color":"#FF3B00"},
    {"name":"7.7x48", "color":"#FF3B00"},
    {"name":"7.7x49", "color":"#FF3B00"},
    {"name":"7.7x50", "color":"#FF3B00"},
    {"name":"7.7x51", "color":"#FF3B00"},
    {"name":"7.7x52", "color":"#FF3B00"},
    {"name":"7.7x53", "color":"#FF3B00"},
    {"name":"7.7x54", "color":"#FF3B00"},
    {"name":"7.7x55", "color":"#FF3B00"},
    {"name":"7.7x56", "color":"#FF3B00"},
    {"name":"7.7x57", "color":"#FF3B00"},
    {"name":"8.8x48", "color":"#FF0000"},
    {"name":"8.8x49", "color":"#FF0000"},
    {"name":"8.8x50", "color":"#FF0000"},
    {"name":"8.8x51", "color":"#FF0000"},
    {"name":"8.8x52", "color":"#FF0000"},
    {"name":"8.8x53", "color":"#FF0000"},
    {"name":"8.8x54", "color":"#FF0000"}
  ]



partial_input =  {"tracks":[
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