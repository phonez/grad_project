## input
	The input is a batch of mol files or mol2 files (deprecated). See phase_02_validation/benchmark/get_ncaa_from_rcsb.ipynb, and there is a script that converts mol object (from smile/inchi, etc) to dipeptide mol file (though few cases get wrong results...). 

## output
	The final output for each ncaa is a params file and a rotlib file.

## scripts
	Run scripts. There are instuctions in each script.

==========
	
## issues
	* Use GetChiralTag() to define L_AA and D_AA...?
	* Inform the version of program we used here...?
	* Name new ncaa with three characters so that it will not conflict with ncaa in rosetta and in my database...