'''
Since MakeRotLib consumes a lot of computing resources, we prepare bash scripts here for run in a cluster.
This script is authored by wfz.
'''

make_rot_lib_app = '/lustre1/chuwang_pkuhpc/rosetta/rosetta_src_2019.47.61047_bundle/main/source/bin/MakeRotLib.mpi.linuxiccrelease'

import os, sys

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

# Each of these bash scripts contains ".lis" with specific phi angle
# Each of lis file contains 36 commands with the same phi angle and different psi angles ranging from -180 to 180.
def gen_hpc_sh(ncaa):
    phi_range = range(-180, 180, 10)
    for phi in phi_range:
        in_list = os.popen(f'ls {ncaa}_rot_lib_options_{phi}_*.in')

        jobs_file = open(f'{ncaa}_rot_lib_{phi}_job.lis','w+')
        hpc_sh_file = open(f'{ncaa}_rot_lib_{phi}_job.sh','w+')

        for in_file in in_list:
            in_file = in_file.replace('\n', '')
            ### -mute MakeRotLibMover and abandon log file (as the number of chi increases, the size of log file will be pretty large) ###
            sh = f'{make_rot_lib_app} -options_file {in_file} -output_logging false -mute protocols.make_rot_lib.MakeRotLibMover -extra_res_fa {ncaa}.params'+'\n'
            jobs_file.write(sh)

        jobs_cmd = f'cat {ncaa}_rot_lib_{phi}_job.lis | parallel -j 20'
        hpc_file.writelines(hpc_sh(f'{ncaa}_rotlib', 1, 20, f'{ncaa}_{phi}', f'{ncaa}_{phi}', jobs_cmd))

        jobs_file.close()
        hpc_sh_file.close()

def main(argv):
    assert[len(argv) == 1]
    ncaa = argv[0]
    gen_hpc_sh(ncaa)

if __name__ == "__main__":
    main(sys.argv[1:])