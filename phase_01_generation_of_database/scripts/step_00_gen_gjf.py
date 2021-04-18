'''
Convert mol2 file (initial structure) to gjf file (gaussian input file).
'''

import sys
import os

in_mol_path = os.path.dirname(sys.path[0]) + "/input/mol/"
gjf_path = os.path.dirname(sys.path[0]) + "/output/gjf/"
chk_path = os.path.dirname(sys.path[0]) + "/output/chk/"
# chk_path = "/home/rotations/zhangf/gen_ncaa_params/output/chk/"

from rdkit import Chem

def get_total_charge(mol):
    total_charge = 0
    for atom in mol.GetAtoms():
        total_charge += atom.GetFormalCharge()
    return total_charge

def get_backbone_atoms(mol):
    backbone = Chem.MolFromSmiles("CC(=O)NCC(=O)NC")  # backbone atoms of dipeptide displayed as ordered
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

def get_phi_atoms_from_backbone(backbone_atoms):
    return (backbone_atoms[5], backbone_atoms[4], backbone_atoms[3], backbone_atoms[1]) # index starts from 0

def get_psi_atoms_from_backbone(backbone_atoms):
    return (backbone_atoms[3], backbone_atoms[4], backbone_atoms[5], backbone_atoms[7])

def gen_gjf_file(gjf_path, chk_path, gjf_name, mol):
    
    gjf_file = gjf_path + gjf_name + ".gjf"
    chk_file = chk_path + gjf_name + ".chk"
    
    with open(gjf_file, "w") as f:
        f.write("%%Chk=%s\n" % (chk_file))
        f.write("#B3LYP/6-31G(d) em=GD3BJ Opt=ModRedundant\n")
        f.write("\n")
        f.write("scan rotamers\n")
        f.write("\n")
        f.write("%d %d\n" %(get_total_charge(mol), 1)) # spin multiplicity is defaulted as 1
        for atom in mol.GetAtoms():
            x = mol.GetConformer().GetAtomPosition(atom.GetIdx()).x
            y = mol.GetConformer().GetAtomPosition(atom.GetIdx()).y
            z = mol.GetConformer().GetAtomPosition(atom.GetIdx()).z
            f.write("%s    %f  %f  %f\n" % (atom.GetSymbol(), x, y, z))
        '''
        # for restricted opitimization: freeze phi and psi
        # however, this syntax is not supported after g09D01
        
        f.write("\n")
        phi_atoms = get_phi_atoms_from_backbone(get_backbone_atoms(mol))
        psi_atoms = get_psi_atoms_from_backbone(get_backbone_atoms(mol))
        f.write("%d %d %d %d  -150.00 F\n" % (phi_atoms[0] + 1, phi_atoms[1] + 1, phi_atoms[2] + 1, phi_atoms[3] + 1))
        f.write("%d %d %d %d  150.00 F\n" % (psi_atoms[0] + 1, psi_atoms[1] + 1, psi_atoms[2] + 1, psi_atoms[3] + 1))
        '''
        f.write("\n")
        f.write("\n")

os.chdir(in_mol_path)
in_mol_list = []

for root, dirs, files in os.walk(in_mol_path):
    for in_mol_file in files:
        if os.path.splitext(in_mol_file)[1] == ".mol": # we should name the ncaa with three characters since now!?
            in_mol_list.append(in_mol_file)

for in_mol_file in in_mol_list:
    gjf_name = os.path.splitext(in_mol_file)[0]
    mol = Chem.MolFromMolFile(in_mol_file, removeHs = False)
    gen_gjf_file(gjf_path, chk_path, gjf_name, mol)