'''
structure optimization (FastRelax) of protein containing ncaa
'''

from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax

ncaa = "BP5"
pdb = "4iww_A"

params_filedir = "/home/zhangf/grad_project/pipeline/params/"
modified_pdb_filedir = "/home/zhangf/grad_project/pipeline/modified_pdb/"
pdb_filedir = "/home/zhangf/grad_project/pipeline/pdb/"

params_file = params_filedir + ncaa + ".params"
rosetta_option = "-in:file:extra_res_fa " + params_file
pdb_file = pdb_filedir + pdb + ".pdb"
modified_pdb_file = modified_pdb_filedir + pdb + ".pdb"

init(extra_options=rosetta_option)
pose = pose_from_file(modified_pdb_file)
origin_pose = pose_from_file(pdb_file)

sfxn = get_score_function(True)
fr = FastRelax(sfxn)
fr.apply(pose)
fr.apply(origin_pose)

print("Before integration:")
print(sfxn(origin_pose))
print("Integrated with ncaa:")
print(sfxn(pose))