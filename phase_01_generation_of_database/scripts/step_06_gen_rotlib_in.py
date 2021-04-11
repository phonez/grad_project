'''
Generate input files (.in) for MakeRotLib application. (to be updated...)
This script is authored by wfz.
'''

def gen_rotlib_in(ncaa, num_chi, sp2_chi):

    phi_range = psi_range = range(-180,180,10)

    for phi in phi_range:
        for psi in psi_range:  

            line_phi = f"PHI_RANGE {phi:>3} {phi:>3} 1"
            line_psi = f"PSI_RANGE {psi:>3} {psi:>3} 1"

            in_name = f"{ncaa}_rot_lib_options_{phi}_{psi}.in"
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
    assert[len(argv) == 3]
    ncaa = argv[0]
    num_chi = (int)(argv[1])
    sp2_chi = (list)(argv[2]) # how to get sp2_chi and num_chi from params.file automatically...?
    gen_rotlib_in(ncaa, num_chi, sp2_chi)

if __name__ == "__main__":
    main(sys.argv[1:])