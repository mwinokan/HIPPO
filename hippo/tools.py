
import re
import numpy as np
from molparse.rdkit import mol_from_smiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import MolFromSmiles, MolToSmiles, AddHs, RemoveHs
import mcol
import mout

from mlog import setup_logger
logger = setup_logger('HIPPO')

from ase.data import chemical_symbols
chemical_symbols = chemical_symbols[1:]
chemical_symbols = sorted(chemical_symbols, key=lambda x: len(x), reverse=True)

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

def flat_inchikey(smiles):
	smiles = sanitise_smiles(smiles)
	return inchikey_from_smiles(smiles)

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

def sanitise_smiles(s, verbosity=False, sanitisation_failed='error', radical='error'):

	assert isinstance(s, str), f'non-string smiles={s}'

	orig_smiles = s

	# if multiple molecules take the largest
	if '.' in s:
		s = sorted(s.split('.'), key=lambda x: len(x))[-1]
	
	# flatten the smiles
	stereo_smiles = s
	smiles = s.replace('@','')
	smiles = smiles.replace('/','')
	smiles = smiles.replace('\\','')

	# remove isotopic stuff
	if smiles_has_isotope(smiles):
		mout.warning(f'Isotope(s) in SMILES: {smiles}')
		smiles = remove_isotopes_from_smiles(smiles)

	# canonicalise
	mol = MolFromSmiles(smiles)
	if mol:
		smiles = MolToSmiles(mol,True)
	elif sanitisation_failed =='error':
		raise SanitisationError
	elif sanitisation_failed =='warning':
		logger.warning(f'sanitisation failed for {smiles=}')

	# check radicals
	reconstruct = False
	for atom in mol.GetAtoms():
		if not atom.GetNumRadicalElectrons():
			continue

		if radical == 'warning':
			logger.warning(f'Radical atom in {smiles=}')
		elif radical == 'error':
			raise SanitisationError(f'Radical atom in {smiles=}')
		elif radical == 'remove':
			logger.warning(f'Removed radical atom')
			atom.SetNumRadicalElectrons(0)         
			smiles = MolToSmiles(mol,True)
			reconstruct = True
			# atom.SetFormalCharge(0)
		else:
			raise NotImplementedError(f'Unknown option {radical=}')

	if reconstruct:
		mol = AddHs(mol)
		mol = RemoveHs(mol, implicitOnly=True)
		smiles = MolToSmiles(mol,True)
		logger.warning(f'New {smiles=}')

	if verbosity:

		if smiles != orig_smiles:

			annotated_smiles_str = orig_smiles.replace('.',f'{mcol.error}{mcol.underline}.{mcol.clear}{mcol.warning}')
			annotated_smiles_str = annotated_smiles_str.replace('@',f'{mcol.error}{mcol.underline}@{mcol.clear}{mcol.warning}')
			
			mout.warning(f'SMILES was changed: {annotated_smiles_str} --> {smiles}')

	return smiles

def sanitise_mol(m):
	from rdkit.Chem import MolToMolBlock, MolFromMolBlock
	return MolFromMolBlock(MolToMolBlock(m))

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

def formula_to_atomtype_dict(formula):

	import re

	remaining = formula.replace('+','')

	atomtype_dict = {}
	
	while remaining:
		for symbol in chemical_symbols:
			if remaining.startswith(symbol):
				remaining = remaining.removeprefix(symbol)

				count_str = re.findall(r'^([0-9]*)', remaining)[0]

				if count_str:
					remaining = remaining.removeprefix(count_str)
					count = int(count_str)
				else:
					count = 1

				atomtype_dict[symbol] = count

				break
		else:
			raise Exception(f'Could not match any chemical symbol to the start of {remaining} {formula=}')

	return atomtype_dict

def combine_atomtype_dicts(atomtype_dicts):

	merged = {}
	for atomtype_dict in atomtype_dicts:
		for symbol, count in atomtype_dict.items():
			if symbol not in merged:
				merged[symbol] = 0

			merged[symbol] += count

	return merged

def atomtype_dict_to_formula(atomtype_dict):
	
	symbols = []

	"""For organic compounds, the order is carbon, hydrogen, 
	then all other elements in alphabetical order of their chemical symbols."""

	key = 'C'
	if key in atomtype_dict:
		value = atomtype_dict[key]
		if value > 1:
			symbols.append(f'C{value}')
		else:
			symbols.append('C')

	key = 'H'
	if key in atomtype_dict:
		value = atomtype_dict[key]
		if value > 1:
			symbols.append(f'{key}{value}')
		else:
			symbols.append(key)

	keys = sorted(list(atomtype_dict.keys()))

	for key in keys:
		value = atomtype_dict[key]
		if key in ['C', 'H']:
			continue
		if value > 1:
			symbols.append(f'{key}{value}')
		else:
			symbols.append(key)

	return ''.join(symbols)

class SanitisationError(Exception):
	...
