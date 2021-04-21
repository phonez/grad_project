## input
	The input is a batch of mol files or mol2 files (deprecated). See phase_02_validation/benchmark/get_ncaa_from_rcsb.ipynb, and there is a script that converts mol object (from smile/inchi, etc) to dipeptide mol file (though few cases get wrong results...). 

## output
	

## scripts
	Run scripts. There are instuctions in each script.

==========
	
## issues
	* Before run hpc_bash, change directory to /output/rotlib (rotlib files are generated in working directory)
	* Use GetChiralTag() to define L_AA and D_AA...?
	* Copy params and rotlib to ncaa_database...
	* Inform the version of program we used here...?
	* For each script, it would make more sense to accomplish the format conversion of a single file (not a batch of files)...