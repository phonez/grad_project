'''
Since MakeRotLib consumes a lot of computing resources, we prepare bash scripts here for run in a cluster.
This script is authored by wfz.
Before run hpc_bash, change directory to /output/rotlib (rotlib files are generated in working directory).
'''

make_rot_lib_app = '/lustre1/chuwang_pkuhpc/rosetta/rosetta_src_2019.47.61047_bundle/main/source/bin/MakeRotLib.mpi.linuxiccrelease'

import os, sys

params_path = os.path.dirname(sys.path[0]) + "/output/params/"
in_path = os.path.dirname(sys.path[0]) + "/output/in/" 
bash_path = os.path.dirname(sys.path[0]) + "/output/bash/" # hpc_bash files

def hpc_sh(J, N, node, o, e, cmd):
    assert type(cmd) == str
    lines = []
    lines.append('#!/bin/bash\n') 
    lines.append(f'#SBATCH -J {J}\n') # Job name
    lines.append('#SBATCH -p cn-long\n')
    lines.append(f'#SBATCH -N {N}\n') # Node number
    lines.append(f'#SBATCH --ntasks-per-node={node}\n')
    lines.append(f'#SBATCH -o {o}.std\n')
    lines.append(f'#SBATCH -e {e}.err\n')
    lines.append('#SBATCH --no-requeue\n')
    lines.append('#SBATCH -A chuwang_g1\n')
    lines.append('#SBATCH --qos=chuwangcnl\n')
    lines.append('\n')
    lines.append('. /apps/source/intel-2018.sh\n')
    lines.append(cmd)
    return lines

# Each of lis file contains 36 commands with the same phi angle and different psi angles ranging from -180 to 180.
def gen_hpc_sh(ncaa):
    phi_range = range(-180, 180, 10)

    hpc_sh_file = open(f'{bash_path}{ncaa}_rot_lib_job.sh','w+')
    jobs_cmd = f'cat {bash_path}{ncaa}_rot_lib_job.lis | parallel -j 20'
    hpc_sh_file.writelines(hpc_sh(f'{ncaa}_rotlib', 1, 20, f'{ncaa}', f'{ncaa}', jobs_cmd))
    hpc_sh_file.close()
    
    jobs_file = open(f'{bash_path}{ncaa}_rot_lib_job.lis','w+')
    for phi in phi_range:
        in_list = os.popen(f'ls {in_path}{ncaa}_rot_lib_options_{phi}_*.in')

        for in_file in in_list:
            in_file = in_file.replace('\n', '')
            ### -mute MakeRotLibMover and abandon log file (as the number of chi increases, the size of log file will be pretty large) ###
            sh = f'{make_rot_lib_app} -options_file {in_file} -output_logging false -mute protocols.make_rot_lib.MakeRotLibMover -extra_res_fa {params_path}{ncaa}.params -use_terminal_residues'+'\n'
            jobs_file.write(sh)
    jobs_file.close()

def main(argv):
    assert[len(argv) == 1]
    ncaa = argv[0]
    gen_hpc_sh(ncaa)

if __name__ == "__main__":
    main(sys.argv[1:])