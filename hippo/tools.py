
import re
import numpy as np
from molparse.rdkit import mol_from_smiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import MolFromSmiles, MolToSmiles
import mcol
import mout

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

def remove_other_ligands(sys, residue_number):
	ligand_residues = [r.number for r in sys['rLIG'] if r.number != residue_number]
	sys.remove_residues_by_index(ligand_residues, fix_indices=False)
	return sys

def inchikey_from_smiles(smiles):
	mol = mol_from_smiles(smiles)
	return MolToInchiKey(mol)

def remove_isotopes_from_smiles(smiles):

	mol = MolFromSmiles(smiles)

	atom_data = [(atom, atom.GetIsotope()) for atom in mol.GetAtoms()]

	for atom, isotope in atom_data:
		if isotope:
			atom.SetIsotope(0)

	return MolToSmiles(mol)

def smiles_has_isotope(smiles, regex=True):
	if regex:
		return re.search(r'([\[][0-9]+[A-Z]+\])',smiles)
	else:
		mol = MolFromSmiles(smiles)
		return any(atom.GetIsotope() for atom in mol.GetAtoms())

def sanitise_smiles(s, verbosity=False):

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

	# canonicalise
	smiles = MolToSmiles(MolFromSmiles(smiles),True)

	if verbosity:

		if smiles != orig_smiles:

			annotated_smiles_str = orig_smiles.replace('.',f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}')
			annotated_smiles_str = annotated_smiles_str.replace('@',f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}')
			
			mout.warning(f'SMILES was changed: {annotated_smiles_str} --> {smiles}')

	return smiles
