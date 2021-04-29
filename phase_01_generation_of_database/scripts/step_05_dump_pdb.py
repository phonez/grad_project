'''
Dump pdb file from params file.
'''

import sys
import os

params_path = os.path.dirname(sys.path[0]) + "/output/params/"
pdb_path = os.path.dirname(sys.path[0]) + "/output/pdb/"

from pyrosetta import *

if len(sys.argv) != 2:
    raise IOError("Input ncaa name with three characters for test!")
ncaa = sys.argv[1]

params_file = params_path + ncaa + ".params"
pdb_file = pdb_path + ncaa + ".pdb"
rosetta_option = "-in:file:extra_res_fa " + params_file

init(extra_options=rosetta_option)
pose_from_sequence(f"X[{ncaa}]").dump_pdb(pdb_file)
