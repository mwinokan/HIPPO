#!/usr/bin/env python3

import os
from db import Database
import molparse as mp
from rdkit import Chem

def main():

	os.system('rm -v test.db')

	mol = mp.rdkit.mol_from_smiles('CCC')
	
	db = Database('test.db')

	db._execute("INSERT INTO COMPOUND(name, smiles, orig_smiles, stereo_smiles, mol) "
		"VALUES(?1, ?2, ?3, ?4, mol_from_smiles(?2))", ('test','c1ccnnc1','c1ccnnc1','c1ccnnc1'))

	# db._execute("INSERT INTO COMPOUND(name, smiles, orig_smiles, stereo_smiles, mol) "
	# 	"VALUES(?1, ?2, ?3, ?4, mol_from_binary_mol(?5))", ('yoyo','CCC','CCC','CCC',mol))

	db._execute("INSERT INTO COMPOUND(name, smiles, orig_smiles, stereo_smiles, mol) "
		"VALUES(?1, ?2, ?3, ?4, mol_from_binary_mol(?5))", ('yoyo','CCC','CCC','CCC',mol.ToBinary()))

	db._execute("SELECT mol_to_binary_mol(mol) FROM COMPOUND")
	mols = [b[0] for b in db.cursor.fetchall()]
	print(mols)
	print([Chem.Mol(b) for b in mols if b is not None])

	# db._execute("INSERT INTO COMPOUND(name, smiles, orig_smiles, stereo_smiles, mol) "
	# 	"VALUES(?1, ?2, ?3, ?4, mol_from_binary_mol(?5))", ('c1cc','c1ccnnc1','c1ccnnc1','c1ccnnc1',mp.rdkit.mol_from_smiles('c1ccnnc1').ToBinary()))

	# db._execute("INSERT INTO COMPOUND(name, smiles, orig_smiles, stereo_smiles, mol) "
	# 	"VALUES(?1, ?2, ?3, ?4, mol_to_binary_mol(?5))", ('waza','CCC','CCC','CCC',mol))

	# db.register_compound(
	# 	name='A', 
	# 	smiles='smiles A', 
	# 	orig_smiles='orig_smiles A', 
	# 	stereo_smiles='stereo_smiles A', 
	# 	mol=mol,
	# )
	# db.register_compound(
	# 	name='B', 
	# 	smiles='smiles B', 
	# 	orig_smiles='orig_smiles B', 
	# 	stereo_smiles='stereo_smiles B', 
	# 	mol=mol,
	# )
	# db.register_compound(
	# 	name='C', 
	# 	smiles='smiles C', 
	# 	orig_smiles='orig_smiles C', 
	# 	stereo_smiles='stereo_smiles C', 
	# 	mol=mol,
	# )

	db._print_table('COMPOUND')

	db.close()

if __name__ == '__main__':
	main()
