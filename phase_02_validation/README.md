## benchmark
	

## stage_01_import_pdb
	This step is to test if a params file has the correct file format. We need to test if rosetta can import a pdb file containing ncaa (with params file from ../ncaa_database).

## stage_02_mutation
	This step is to test if a residue at specified position in protein can be mutated to a ncaa.
	
## stage_03_optimization
	This step is to test if a ncaa residue in protein has a reasonable conformation after structure optimization as well as a params file with correct bonds and internal coordinates.
	
## stage_04_regeneration
	This step is to test if our protocol can regenerate