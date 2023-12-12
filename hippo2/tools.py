
import mout
import numpy as np
from rdkit import Chem
import re

def df_row_to_dict(df_row):

	assert len(df_row) == 1

	data = {}

	for col in df_row.columns:

		if col == 'Unnamed: 0':
			continue

		value = df_row[col].values[0]

		if not isinstance(value,str) and np.isnan(value):
			value = None

		data[col] = value

	return data

def remove_isotopes_from_smiles(smiles):

	mol = Chem.MolFromSmiles(smiles)

	atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]

	for atom, isotope in atom_data:
		if isotope:
			atom.SetIsotope(0)

	return Chem.MolToSmiles(mol)

def smiles_has_isotope(smiles, regex=True):
	if regex:
		return re.search(r'([\[][0-9]+[A-Z]+\])',smiles)
	else:
		mol = Chem.MolFromSmiles(smiles)
		return any(atom.GetIsotope() for atom in mol.GetAtoms())
