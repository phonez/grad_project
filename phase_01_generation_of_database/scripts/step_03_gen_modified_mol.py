'''
Convert mol file to modified mol file. (to be updated...)
'''

import sys
import os

mol_path = os.path.dirname(sys.path[0]) + "/output/mol/"
modified_mol_path = os.path.dirname(sys.path[0]) +"/output/modified_mol/"

from rdkit import Chem

def get_total_charge(mol):
    total_charge = 0
    for atom in mol.GetAtoms():
        total_charge += atom.GetFormalCharge()
    return total_charge

def get_backbone_atoms(mol):
    backbone = Chem.MolFromSmiles("CC(=O)NCC(=O)NC")  # backbone atoms of mol displayed as ordered
    backbone_atoms = mol.GetSubstructMatch(backbone)
    num_Hs_in_ACE = 0
    num_Hs_in_NME = 0
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[0]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom covalent with C atom in ACE 
            num_Hs_in_ACE += 1
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[8]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom covalent with C atom in NME 
            num_Hs_in_NME += 1
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[7]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom covalent with N atom in NME
            num_Hs_in_NME += 1
    assert num_Hs_in_ACE == 3 and num_Hs_in_NME == 4
    return backbone_atoms

def get_carbon_in_ACE(backbone_atoms): # POLY_LOWER ATOM in modified mol file
    return backbone_atoms[1] + 1

def get_nitrogen_in_NME(backbone_atoms): # POLY_UPPER ATOM in modified mol file
    return backbone_atoms[7] + 1

def get_capping_atoms(mol): # except POLY_LOWER ATOM and POLY_UPPER ATOM in ACE and NME
    backbone_atoms = get_backbone_atoms(mol)
    capping_atoms = set()
    capping_atoms.add(backbone_atoms[0] + 1)
    capping_atoms.add(backbone_atoms[8] + 1)   
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[0]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom
            capping_atoms.add(neighbor.GetIdx() + 1)
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[8]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom
            capping_atoms.add(neighbor.GetIdx() + 1)
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[1]).GetNeighbors():
        if neighbor.GetAtomicNum() == 8: # oxygen atom
            capping_atoms.add(neighbor.GetIdx() + 1)
    for neighbor in mol.GetAtomWithIdx(backbone_atoms[7]).GetNeighbors():
        if neighbor.GetAtomicNum() == 1: # hydrogen atom
            capping_atoms.add(neighbor.GetIdx() + 1)
    return capping_atoms

def is_aromatic(mol):
    pass

os.chdir(mol_path)
mol_list = []

for root, dirs, files in os.walk(mol_path):
    for mol_file in files:
        if os.path.splitext(mol_file)[1] == ".mol":
            mol_list.append(mol_file)

for mol_file in mol_list:
    mol_name = os.path.splitext(mol_file)[0]  
    mol_file = mol_path + mol_name + ".mol"
    modified_mol_file = modified_mol_path + mol_name + ".mol"
    
    mol = Chem.MolFromMolFile(mol_file, removeHs = False)
    backbone_atoms = get_backbone_atoms(mol)
    
    ROOT = backbone_atoms[3] + 1
    POLY_N_BB = backbone_atoms[3] + 1
    POLY_CA_BB = backbone_atoms[4] + 1
    POLY_C_BB = backbone_atoms[5] + 1
    POLY_O_BB = backbone_atoms[6] + 1
    POLY_IGNORE = get_capping_atoms(mol)
    POLY_UPPER = get_nitrogen_in_NME(backbone_atoms)
    POLY_LOWER = get_carbon_in_ACE(backbone_atoms)
    POLY_CHARGE = get_total_charge(mol)
    POLY_PROPERTIES = ["PROTEIN"] # to be updated...
    
    with open(mol_file, "r") as f_read, open(modified_mol_file, "w") as f:
        read_lines = f_read.readlines()
        cursor = -1
        while (read_lines[cursor] == "\n" or read_lines[cursor][0] == 'M'):
            cursor = cursor - 1
        cursor = cursor + 1

        # we need to adjust the relative position of records of bond info
        # to span atom tree correctly in molfile_to_params_polymer.py (due to some stupid functions...)
        # N -> CA -> C -> UPPER
        # and we'll avoid the confusing chi bug "N CA C O"
        num_atoms = int(read_lines[3][0:3]) # refers to scripts/python/rosetta_py
        first_bond_idx = 3 + num_atoms + 1

        N_CA_line = ""
        CA_C_line = ""
        C_UPPER_line = ""
        for line in read_lines:
            if len(line.split()) > 1:
                tmp = {line.split()[0], line.split()[1]} 
                if {f"{POLY_N_BB}", f"{POLY_CA_BB}"} == tmp:
                    N_CA_line = line
                if {f"{POLY_CA_BB}", f"{POLY_C_BB}"} == tmp:
                    CA_C_line = line
                if {f"{POLY_C_BB}", f"{POLY_UPPER}"} == tmp:
                    C_UPPER_line = line
        
        # swap position
        read_lines.remove(N_CA_line)
        read_lines.remove(CA_C_line)
        read_lines.remove(C_UPPER_line)
        read_lines.insert(first_bond_idx, C_UPPER_line)
        read_lines.insert(first_bond_idx, CA_C_line)
        read_lines.insert(first_bond_idx, N_CA_line)
        
        f.writelines(read_lines[:cursor])
        f.write("M  ROOT %d\n" % (ROOT))
        f.write("M  POLY_N_BB %d\n" % (POLY_N_BB))
        f.write("M  POLY_CA_BB %d\n" % (POLY_CA_BB))
        f.write("M  POLY_C_BB %d\n" % (POLY_C_BB))
        f.write("M  POLY_O_BB %d\n" % (POLY_O_BB))
        f.write("M  POLY_IGNORE")
        for i in POLY_IGNORE:
            f.write(" %d" % (i))
        f.write("\n")
        f.write("M  POLY_UPPER %d\n" % (POLY_UPPER))
        f.write("M  POLY_LOWER %d\n" % (POLY_LOWER))
        f.write("M  POLY_CHG %d\n" % (POLY_CHARGE))
        f.write("M  POLY_PROPERTIES")
        for str in POLY_PROPERTIES:
            f.write(" %s" % (str))
        f.write("\n")
        f.write("M  END\n")