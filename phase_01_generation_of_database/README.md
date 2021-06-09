## input
	The input is a batch of mol files or mol2 files (deprecated). 
├─input
   ├─mol
   └─mol2

## output
	The final output for each ncaa is a params file and a rotlib file.
├─output
   ├─bash
   ├─chk
   ├─gjf
   ├─in
   ├─log
   ├─modified_mol
   ├─mol
   ├─params
   ├─pdb
   └─rotlib

## scripts
	Run scripts. There are instuctions in each script.


	
## issues
	* Use GetChiralTag() to define L_AA and D_AA...?
	* Name new ncaa with three characters so that it will not conflict with ncaa in rosetta and in my database...