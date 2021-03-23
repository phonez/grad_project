'''
pdb file --> modified pdb file
'''

params_filedir = "/home/zhangf/grad_project/pipeline/params/"
pdb_filedir = "/home/zhangf/grad_project/pipeline/pdb/"
modified_pdb_filedir = "/home/zhangf/grad_project/pipeline/modified_pdb/"

from pyrosetta import *
from rdkit import Chem

# generate ncaa.pdb file from params file
init(extra_options="-in:file:extra_res_fa /home/zhangf/grad_project/sampling-score_sys/params/BP5.params")
BP5_pose=pose_from_sequence("X[BP5]")
BP5_pose.dump_pdb("BP5_from_rosetta.pdb")
BP5_from_rosetta = Chem.MolFromPDBFile("../convert_params_to_pdb/BP5_from_rosetta.pdb")

# generate ncaa.pdb file from pdb file (with water, sulfate, and another identical chain all removed previously)
get_BP5_from_pdb = [] 
with open('./4iww_A.pdb', "r") as f_read_pdb, open('./BP5_from_pdb.pdb', "w") as f_write:
    lines_in_f_read_pdb = f_read_pdb.readlines()
    for line in lines_in_f_read_pdb:
        if line[0:6] == "HETATM" and line[17:20] == "BP5": # 找到pdb文件中记录ncaa的文本（可能有多个ncaa）
            BP5_from_pdb.append(line)
    f_write.writelines(get_BP5_from_pdb)
BP5_from_pdb = Chem.MolFromPDBFile("./BP5_from_pdb.pdb") # 从含有ncaa的蛋白质截取下来的ncaa

# get corresponding relationship of atoms between ncaa_from_rosetta and ncaa_from_pdb
BP5 = BP5_from_rosetta.GetSubstructMatch(BP5_from_pdb)
dict = {}
for i in range(BP5_from_rosetta.GetNumAtoms()):
    dict[i] = BP5[i]

# modify pdb file
atom_symbol_from_rosetta = []
with open('../convert_params_to_pdb/BP5_from_rosetta.pdb', "r") as f_read_rosetta, open('./4iww_A.pdb', "r") as f_read_pdb, open('./modified_pdb', "w") as f_write:
        lines_in_f_read_rosetta = f_read_rosetta.readlines()
        for line in lines_in_f_read_rosetta:
            if line[0:6] == "HETATM" and line[17:20] == "BP5" and line[13] !="H":
                atom_symbol_from_rosetta.append(line[12:16])

        lines_in_f_read_pdb = f_read_pdb.readlines()
        modified_pdb_lines = []
        ncaa_lines_num = 0
        for line in lines_in_f_read_pdb:
            if ncaa_lines_num == BP5_from_rosetta.GetNumAtoms():
                ncaa_lines_num = 0
            
            if line[0:6] == "HETATM" and line[17:20] == "BP5": # 找到pdb文件中记录ncaa的文本（可能有多个ncaa）
                modified_line = "ATOM  " + line[6:12] + atom_symbol_from_rosetta[dict[ncaa_lines_num]] + line[16:]
                modified_pdb_lines.append(modified_line)
                ncaa_lines_num = ncaa_lines_num + 1
            else:
                modified_pdb_lines.append(line)

        f_write.writelines(modified_pdb_lines)