## benchmark
	Select proteins or protein-peptide complexes containing ncaa from pdb as our benchmarks. They are tested from stage_01 to stage_03.

## stage_01_import_pdb
	This step is to test if rosetta can import the pdb file containing a specific ncaa (with params file from ../ncaa_database).

## stage_02_mutation
	This step is to test if the residue in protein can be mutated to a specific ncaa.
	
## stage_03_optimization
	This step is to test if a specific ncaa residue in protein has the similar conformation compared to native structure after mutation and structure optimization.
	
## stage_04_regeneration
	This step is to test if our protocol can regenerate results consistent with experimental results (binding affinity of protein-peptide complex) by point mutation scan.