## Dependencies

These scripts are tested to work under Python 3.6, and Python 2.7 should also be installed (with command python2 in bash env).

Required dependencies are 
-   [rdkit](https://rdkit.org/), version=2020.09.1.

Optional dependencies are 
-   [pyrosetta](https://www.pyrosetta.org/), version=2021.11+release.e9f4797
-   [Gaussian](http://gaussian.com/), version=g09D01
-   [openbabel](http://openbabel.org/), version=3.1.0

## Auto-generation of params file from SMILES

### Input
* See example.txt. The first column represents the 3-character name, the second column represents the SMILES txt.
* A directory containing mol files with given dipeptide 3D structure.

### Generate params file
To generate params file, use gen_all.py script with the following options:
* `input`             str, name of a file storing SMILES of NCAA **or** mol files directory. See [Input](###Input).
* `gaussian_opt`      bool, if True, use Gaussian 09 for the structure optimization of a dipeptide molecule; if False, use MMFF force field for optimization in rdkit from ./scripts/step_00_gen_dipeptide.py. False as default.
* `check_params`      bool, if True, use pyrosetta to dump pdb to check if params file works. False as default.
* `makerotlib`        str, name of the path of MakeRotLib app for the generation of rotlib file. An empty string as default.
* `clean`             bool, if True, remove intermediate files. True as default.

```
cd scripts
python gen_all.py --input=../example.txt --clean
cd ../output/params
```
Or if you already have the dipeptide mol files in ./mol (See [Note](###Note)):
```
cd scripts
python gen_all.py --input=../mol/ --clean
cd ../output/params
```
The params file is in the current directory.

### Note
- Since there could be errors for a few NCAAs in ./scripts/step_00_gen_dipeptide.py and ./scripts/step_04_molfile_to_params_polymer.py, checking the structure of dipeptide or dumped pdb is recommended.
    -   Sometimes rdkit (in step_00_gen_dipeptide.py) can't find the 3D conformation with lowest energy from SMILES (eg. kew-boat conformation but not chair conformation), and gaussian optimization may not converge. In this case, we could increase the `numconfs` in AllChem.EmbedMultipleConfs function.
    -   In few cases the structure of dipeptide would be wrong due to the defect of the algorithm of adding NME and ACE group. In this case, one could input a mol file with a given dipeptide structure to skip this gen_dipeptide step.
- Using gaussian for optimization and making rotlib file are both time-consuming, modify the corresponding scripts according to your computing resource, if necessary.
