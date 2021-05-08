'''
Convert pdb file to modified pdb file (which can be imported by rosetta with corresponding params file). (to be updated..)
Before conversion, we should manully clean pdb file , e.g., remove solvent, remove small molecule, leave one chain of poly-
mers (if possible), and remove ANISOU and alternate location indicator records in pdb file.. 
'''

import os, sys

params_path = "/home/zhangf/grad_project/ncaa_database/params/"
pdb_path = os.path.dirname(sys.path[0]) + "/input/"
modified_pdb_path = os.path.dirname(sys.path[0]) + "/output/"

from pyrosetta import *
from rdkit import Chem

# generate ncaa.pdb file from params file
# what if there is ncaa already in rosetta database (to be updated...)
def gen_ncaa_pdb_from_rosetta(ncaa): 
    params_file = params_path + ncaa + ".params"
    rosetta_option = "-in:file:extra_res_fa " + params_file
    init(extra_options=rosetta_option)
    ncaa_pose = pose_from_sequence("X[" + ncaa + "]")
    ncaa_pdb_from_rosetta =  pdb_path + ncaa + "_from_rosetta.pdb"
    ncaa_pose.dump_pdb(ncaa_pdb_from_rosetta)

# generate ncaa.pdb file from pdb file
# the name of ncaa in pdb is consistent with that in rosetta
def gen_ncaa_pdb_from_pdb(pdb, ncaa):
    get_ncaa_from_pdb = []
    ncaa_pdb_from_pdb =  pdb_path + ncaa + "_from_pdb.pdb"
    pdb_file = pdb_path + pdb + ".pdb"

    ncaa_pos = ["X", -1]
    with open(pdb_file, "r") as f_read_pdb, open(ncaa_pdb_from_pdb, "w") as f_write:
        lines_in_f_read_pdb = f_read_pdb.readlines()
        for line in lines_in_f_read_pdb:
            cur_pos = []
            if len(line) >= 80: # for line that records ATOM or HETATM
                cur_pos = [line[21], (int)(line[22:26])] # pdb format, chain name and residue position
            if (ncaa_pos[1] < 0 or ncaa_pos == cur_pos) and line[0:6] == "HETATM" and line[17:20] == ncaa: # avoid duplication
                ncaa_pos = cur_pos
                get_ncaa_from_pdb.append(line[:78] + "  \n") # pdb from rosetta doesn't contain atom charge, so we have to discard it; else a problem will occur in substructure match
        f_write.writelines(get_ncaa_from_pdb)

# get corresponding relationship of atoms between ncaa_from_rosetta and ncaa_from_pdb
def get_atom_dict(ncaa): 
    ncaa_pdb_from_pdb =  pdb_path + ncaa + "_from_pdb.pdb"
    ncaa_pdb_from_rosetta =  pdb_path + ncaa + "_from_rosetta.pdb"
    ncaa_from_pdb = Chem.MolFromPDBFile(ncaa_pdb_from_pdb)
    ncaa_from_rosetta = Chem.MolFromPDBFile(ncaa_pdb_from_rosetta)
    ncaa = ncaa_from_rosetta.GetSubstructMatch(ncaa_from_pdb)
    dict = {}
    for i in range(ncaa_from_rosetta.GetNumAtoms()):
        dict[i] = ncaa[i]
    return dict

# modify pdb file
def modify_pdb_file(pdb, ncaa):
    ncaa_pdb_from_rosetta =  pdb_path + ncaa + "_from_rosetta.pdb"
    ncaa_from_rosetta = Chem.MolFromPDBFile(ncaa_pdb_from_rosetta)
    pdb_file = pdb_path + pdb + ".pdb"
    modified_pdb_file = modified_pdb_path + pdb + ".pdb"
    
    atom_dict = get_atom_dict(ncaa)
    atom_symbol_from_rosetta = [] # CA, CB, ect.
    with open(ncaa_pdb_from_rosetta, "r") as f_read_rosetta, open(pdb_file, "r") as f_read_pdb, open(modified_pdb_file, "w") as f_write:
            lines_in_f_read_rosetta = f_read_rosetta.readlines()
            for line in lines_in_f_read_rosetta:
                if line[0:6] == "HETATM" and line[17:20] == ncaa and line[13] !="H":
                    atom_symbol_from_rosetta.append(line[12:16])

            lines_in_f_read_pdb = f_read_pdb.readlines()
            modified_pdb_lines = []
            ncaa_lines_num = 0
            for index, line in enumerate(lines_in_f_read_pdb):
                if ncaa_lines_num == ncaa_from_rosetta.GetNumAtoms():
                    ncaa_lines_num = 0
                if line[0:3] == "TER": # connet ncaa with canonical aa (avoid unexpected termination)
                    next_line = lines_in_f_read_pdb[index + 1]
                    last_line = lines_in_f_read_pdb[index - 1]
                    if next_line[0:6] == "HETATM" and next_line[17:20] == ncaa and next_line[21] == last_line[21]:
                        continue
                if line[0:6] == "HETATM" and line[17:20] == ncaa: 
                    if line[13] !="H":
                        modified_line = "ATOM  " + line[6:12] + atom_symbol_from_rosetta[atom_dict[ncaa_lines_num]] + line[16:]
                    else:
                        modified_line = "ATOM  " + line[6:]
                    modified_pdb_lines.append(modified_line)
                    ncaa_lines_num = ncaa_lines_num + 1
                else:
                    modified_pdb_lines.append(line)
            f_write.writelines(modified_pdb_lines)

def main(argv):
    assert[len(argv) == 2]
    gen_ncaa_pdb_from_rosetta(argv[1]) # ncaa to be validated
    gen_ncaa_pdb_from_pdb(argv[0], argv[1]) # benchmark pdb
    modify_pdb_file(argv[0], argv[1])

if __name__ == "__main__":
    main(sys.argv[1:])