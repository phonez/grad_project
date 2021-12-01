import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input", type=str, required=True, default="", help="Smiles list file path or mol files directory.")
parser.add_argument("--gaussian_opt", action='store_true', help="Use gaussian 09 for structure optimization of molecule, else we use dipeptides from step_00_gen_dipeptide.py with MMFF force field optimization.")
parser.add_argument("--check_params", action='store_true', help="Dump pdb through pyrosetta to check if params file works.")
parser.add_argument("--makerotlib", type=str, default="", help="Make rotlib or not, with path of MakeRotLib app as input. This is a time-consuming step, usually taken in hpc cluster.")
parser.add_argument("--clean", action='store_true', help="Remove intermediate files, with params file remained.")

args = parser.parse_args()
cmds = []
if os.path.isdir(args.input):
    cmds.append(f"cp {args.input}/*.mol ../input/mol/")
else:
    cmds.append(f"python step_00_gen_dipeptide.py {args.input}")
if args.gaussian_opt:
    cmds.append("python step_01_gen_gjf.py && bash step_01_run_gaussian.sh")
    cmds.append("bash step_02_gen_mol.sh")
else:
    cmds.append("cp ../input/mol/* ../output/mol/")
cmds.append("python step_03_gen_modified_mol.py")
cmds.append("bash step_04_gen_params.sh")
if args.check_params:
    cmds.append("bash step_05_check_params.sh")
if len(args.makerotlib) != 0:
    cmds.append(f"python step_06_gen_rotlib_in.py && python step_06_gen_hpc_bash.py {args.makerotlib}")
if args.clean:
    cmds.append("rm -r ../input/ && cd ../output/ && rm -r gjf chk log mol modified_mol")
    if not args.check_params:
        cmds.append("cd ../output/ && rm -r pdb")
    if len(args.makerotlib) == 0:
        cmds.append("cd ../output/ && rm -r in bash rotlib")

os.makedirs("../input/")
os.makedirs("../input/mol/")
os.makedirs("../output/")
os.makedirs("../output/gjf/")
os.makedirs("../output/chk/")
os.makedirs("../output/log/")
os.makedirs("../output/mol/")
os.makedirs("../output/modified_mol/")
os.makedirs("../output/params/")
os.makedirs("../output/pdb/")
os.makedirs("../output/in/")
os.makedirs("../output/bash/")
os.makedirs("../output/rotlib/")

for cmd in cmds:
    print(f"Current job: {cmd}")
    os.system(cmd)