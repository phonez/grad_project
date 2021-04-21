'''
Generate input files (.in) for MakeRotLib application. (to be updated...)
This script is authored by wfz.
'''

import os
import sys

params_path = os.path.dirname(sys.path[0]) + "/output/params/"
in_path = os.path.dirname(sys.path[0]) + "/output/in/" # input files for generation fo rotlib

def get_chi_info(ncaa):
    log_file = params_path + ncaa + ".log"
    params_file = params_path + ncaa + ".params"
    chi_info = []
    proton_chis = set()
    with open(params_file, "r") as f:
        for line in f.readlines():
            if line[:10] == "PROTON_CHI": # discard proton chi
                proton_chis.add(line[11])
    with open(log_file, "r") as f:
        for line in f.readlines():
            if line[4] not in proton_chis and line[:3] == "CHI": # make sure there isn't wrong chis in params file or log file previously
                (chi_num, chi_type) = ((int)(line[4]), line[-8:].strip())
                chi_info.append((chi_num, chi_type))
    if len(chi_info) > 4:
        chi_info = chi_info[:4] # discard other chis barbarically...rotlib can only deal with no more than 4 chis
    # get the number of chis and the position of sp3-sp2 chi
    num_chis = len(chi_info)
    sp3_sp2_chi_num = []
    for index, chi in enumerate(chi_info):
        if chi[1] == 'sp3-sp2' or chi[1] == 'sp2-sp2': # as for sp2-sp2 bond, sample 6 bins as well (more is better(we are stupid in chemistry, and such case is less frequent)...in fact, for such bond in BP5 and 4AF, 2 bins (30/60) is enough)
            sp3_sp2_chi_num.append(index + 1) # ...assume that MakeRotLib app gets the chis as the order presented in params file
    return (num_chis, sp3_sp2_chi_num)

def gen_rotlib_in(ncaa):
    (num_chi, sp2_chi) = get_chi_info(ncaa)

    phi_range = psi_range = range(-180,180,10)
    for phi in phi_range:
        for psi in psi_range:  

            line_phi = f"PHI_RANGE {phi:>3} {phi:>3} 1"
            line_psi = f"PSI_RANGE {psi:>3} {psi:>3} 1"

            in_name = f"{in_path}{ncaa}_rot_lib_options_{phi}_{psi}.in"
            with open(in_name, 'w+') as f:
                f.write(f'AA_NAME {ncaa}\n')
                f.write('OMG_RANGE 180 180 1\n')
                f.write(line_phi + '\n')
                f.write(line_psi + '\n')
                f.write('EPS_RANGE 180 180 1\n')
                f.write(f'NUM_CHI {num_chi}\n')
                f.write('NUM_BB 2\n')

                for chi in range(1, num_chi + 1):
                    f.write(f'CHI_RANGE {chi} 0  330  30\n')

                for rotwells in range(1, num_chi + 1):
                    if rotwells in sp2_chi:
                        f.write(f'ROTWELLS {rotwells} 6  180 210 240 270 300 330\n')
                    else:
                        f.write(f'ROTWELLS {rotwells} 3  60 180 300\n')

def main(argv):
    assert[len(argv) == 1]
    ncaa = argv[0]
    gen_rotlib_in(ncaa)

if __name__ == "__main__":
    main(sys.argv[1:])