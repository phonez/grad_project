'''
Calculate chi angles of each chi of specified ncaa (probably more than one) in protein.
'''

from Bio.PDB import *
import math

greek_alphabet = ['A', 'B', 'G', 'D', 'E', 'Z', 'T', 'I', 'K', 'L', 'M', 'N', 'X', 'O', 'P', 'R', 'S', 'U', 'P', 'C']

def comp(elem):
    idx = -1
    for index, greek in enumerate(greek_alphabet):
        if greek == elem[1][1]:
            idx = index
            break
    return idx

def get_chi_atoms(ncaa):
    params_file = f"/home/zhangf/grad_project/ncaa_database/params/{ncaa}.params"
    
    chi_atoms = []
    with open(params_file, "r") as f:
        for line in f.readlines():
            if line[:3] == "CHI":
                atoms = line.strip().split()[-4:]
                chi_atoms.append(atoms)
            if line[:10] == "PROTON_CHI": # discard proton chi
                chi_atoms.pop()

    chi_atoms.sort(key=comp) # sorted as order in atom spanning tree
    return chi_atoms

def get_chi_angles(ncaa, pdb_file):
    ''' Get all chi angles in each specified ncaa in protein. '''

    p = PDBParser()
    protein = p.get_structure("X", pdb_file)
    ncaa_res = []
    for res in protein.get_residues():
        if res.get_resname() == ncaa:
            ncaa_res.append(res)
    # print(f"This protein has {len(ncaa_res)} position(s) inserted with {ncaa}.")

    ncaa_chi_angles = []
    for res in ncaa_res:
        chi_angles = []
        for chi_atoms in get_chi_atoms(ncaa):
            [atom_1, atom_2, atom_3, atom_4] = [res[i] for i in chi_atoms]
            angle = calc_dihedral(atom_1.get_vector(), atom_2.get_vector(), atom_3.get_vector(), atom_4.get_vector())
            chi_angles.append(round(angle * 180 / math.pi, 2))
        ncaa_chi_angles.append(chi_angles)
    print(f"After optimization, chi angles in each ncaa are: {ncaa_chi_angles}")

    return ncaa_chi_angles