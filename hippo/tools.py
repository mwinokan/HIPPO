
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

def remove_other_ligands(sys, residue_number, chain):
	
	ligand_residues = [r.number for r in sys['rLIG'] if r.number != residue_number]
	
	# if ligand_residues:
	for c in sys.chains:
		if c.name != chain:
			c.remove_residues(names=['LIG'], verbosity=0)
		elif ligand_residues:
			c.remove_residues(numbers=ligand_residues, verbosity=0)

	# print([r.name_number_str for r in sys['rLIG']])

	assert len([r.name_number_str for r in sys['rLIG']]) == 1, f"{sys.name} {[r.name_number_str for r in sys['rLIG']]}"

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

def pose_gap(a, b):

	from numpy.linalg import norm
	from molparse.rdkit import mol_to_AtomGroup

	min_dist = None

	a = mol_to_AtomGroup(a.mol)
	b = mol_to_AtomGroup(b.mol)

	for atom1 in a.atoms:
		for atom2 in b.atoms:
			dist = norm(atom1.np_pos - atom2.np_pos)
			if min_dist is None or dist < min_dist:
				min_dist = dist

	return min_dist
