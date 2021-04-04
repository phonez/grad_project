'''
initial structure mol2 file --> gjf file 
'''

mol2_filedir = "/home/zhangf/grad_project/pipeline/mol2/"
gjf_filedir = "/home/zhangf/grad_project/pipeline/gjf/"
chk_filedir = "/home/zhangf/grad_project/pipeline/chk/"

from rdkit import Chem
import os

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
    return (backbone_atoms[5], backbone_atoms[4], backbone_atoms[3], backbone_atoms[1]) # index starts from 0!?

def get_psi_atoms_from_backbone(backbone_atoms):
    return (backbone_atoms[3], backbone_atoms[4], backbone_atoms[5], backbone_atoms[7])

def gen_gjf_file(gjf_filedir, chk_filedir, gjf_filename, mol):
    
    gjf_file = gjf_filedir + gjf_filename + ".gjf"
    chk_file = chk_filedir + gjf_filename + ".chk"
    
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
        
        f.write("\n")
        phi_atoms = get_phi_atoms_from_backbone(get_backbone_atoms(mol))
        psi_atoms = get_psi_atoms_from_backbone(get_backbone_atoms(mol))
        f.write("%d %d %d %d  -150.00 F\n" % (phi_atoms[0] + 1, phi_atoms[1] + 1, phi_atoms[2] + 1, phi_atoms[3] + 1))
        f.write("%d %d %d %d  150.00 F\n" % (psi_atoms[0] + 1, psi_atoms[1] + 1, psi_atoms[2] + 1, psi_atoms[3] + 1))
        '''
        f.write("\n")
        f.write("\n")

os.chdir(mol2_filedir)
mol2_filelist = []

for root, dirs, files in os.walk(mol2_filedir):
    for file in files:
        if os.path.splitext(file)[1] == ".mol2":
            mol2_filelist.append(file)

for mol2_file in mol2_filelist:
    gjf_filename = os.path.splitext(mol2_file)[0]
    mol = Chem.MolFromMol2File(mol2_file, removeHs = False)
    gen_gjf_file(gjf_filedir, chk_filedir, gjf_filename, mol)