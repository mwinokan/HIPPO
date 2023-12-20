
import mout
import numpy as np
from rdkit import Chem
import re
import mcol

def df_row_to_dict(df_row):

	assert len(df_row) == 1, f'{len(df_row)=}'

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

def clean_smiles(s, verbosity=0):

	orig_smiles = s

	# if multiple molecules take the largest
	if '.' in s:
		s = sorted(s.split('.'), key=lambda x: len(x))[-1]
	
	# flatten the smiles
	stereo_smiles = s
	smiles = s.replace('@','')

	# remove isotopic stuff
	if smiles_has_isotope(smiles):
		mout.warning(f'Isotope(s) in SMILES: {smiles}')
		smiles = remove_isotopes_from_smiles(smiles)

	if verbosity:

		if smiles != orig_smiles:

			annotated_smiles_str = orig_smiles.replace('.',f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}')
			annotated_smiles_str = annotated_smiles_str.replace('@',f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}')
			
			mout.warning(f'SMILES was changed: {annotated_smiles_str} --> {smiles}')

	# canonicalise
	smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles),True)

	return dict(smiles=smiles, orig_smiles=orig_smiles, stereo_smiles=stereo_smiles)

def number_to_base(n, b):
	if n == 0:
		return [0]
	digits = []
	while n:
		digits.append(int(n%b))
		n //= b
	return digits[::-1]

def alphabet_index(index, size=3):

	if index >= pow(26,size):
		raise Exception('index too big for string length')

	digits = number_to_base(index, 26)

	for i in range(size - len(digits)):
		digits = [0] + digits
	
	return ''.join([chr(65+d) for d in digits])
   