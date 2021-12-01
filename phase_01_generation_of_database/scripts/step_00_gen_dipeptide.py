import sys
from rdkit import Chem
from rdkit.Chem import AllChem
# from rdkit.Chem import Draw

# load ncaa from smiles
ncaa_info = []
smiles_file = sys.argv[1]
with open(smiles_file, "r") as f:
    records = f.readlines()
    for record in records:
        [ncaa_entry, smiles] = record.strip().split()
        ncaa_mol = Chem.MolFromSmiles(smiles)
        ncaa_info.append([ncaa_entry, ncaa_mol])

# add NME and ACE
def modify_mol(ncaa_mol):
    dipep_bb_mol = Chem.MolFromSmarts("CC(=O)NCC(=O)NC")
    nme_mol = Chem.MolFromSmiles("C(=O)NC")
    ace_mol = Chem.MolFromSmiles("NC(=O)C")
    amino_mol = Chem.MolFromSmarts("N")
    acid_mol = Chem.MolFromSmarts("C(=O)O")

    add_nme_mols = AllChem.ReplaceSubstructs(ncaa_mol, acid_mol, nme_mol, useChirality=True)
    for mol_1 in add_nme_mols:
        add_ace_mols = AllChem.ReplaceSubstructs(mol_1, amino_mol, ace_mol, useChirality=True)  
        for mol_2 in add_ace_mols:
            if mol_2.HasSubstructMatch(dipep_bb_mol):
                Chem.SanitizeMol(mol_2)
                return mol_2

modified_ncaa_mols = []
for ncaa_entry, ncaa_mol in ncaa_info:
    mol = modify_mol(ncaa_mol)
    modified_ncaa_mols.append([mol, ncaa_entry])
### check structure, because this script is not compatible with a few cases ###
# image = Draw.MolsToGridImage([x[0] for x in modified_ncaa_mols], molsPerRow=3, subImgSize=(300,300), legends=[x[1] for x in modified_ncaa_mols], returnPNG=False, maxMols=100)

# 2D to 3D
mol_path = "../input/mol/"
for mol, ncaa_entry in modified_ncaa_mols:
    m = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m)
    AllChem.EmbedMultipleConfs(m, numConfs=400)
    AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0)

    m.SetProp("_Name", ncaa_entry)
    ncaa_file = mol_path + ncaa_entry + ".mol"
    print(Chem.MolToMolBlock(m), file=open(ncaa_file, 'w+'))