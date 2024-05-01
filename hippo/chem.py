
from mout import debug_log

import logging
logger = logging.getLogger('HIPPO')

"""

Checks
======

- Num heavy atoms difference
- Formula checks
- Num rings difference

"""

SUPPORTED_CHEMISTRY = {

	"Amidation" : 
	{
		"heavy_atoms_diff":1, 
		"rings_diff":0,
		"atomtype":
		{
			"removed" : { "O":1, "H":2 },
		},
	},

	"Williamson_ether_synthesis" : 
	{
		"heavy_atoms_diff":1, 
		"rings_diff":0,
		"atomtype":
		{
			"removed" : { "Ha":1, "H":1 }, # any halogen
		},
	},

	"N-Boc_deprotection" : 
	{
		"heavy_atoms_diff":7, 
		"rings_diff":0,
		"atomtype":
		{
			"removed" : { "O":2, "C":5, "H":8 },
		},
	},

	"TBS_alcohol_deprotection" : 
	{
		"heavy_atoms_diff":7, 
		"rings_diff":0,
		"atomtype":
		{
			"removed" : { "C":6, "Si":1, "H":14 },
		},
	},

	"Sp3-sp2_Suzuki_coupling" : 
	{
		"heavy_atoms_diff":10, 
		"rings_diff":1,
		"atomtype":
		{
			"removed" : { "C":6, "O":2, "B":1, "Ha":1, "H":12 }, # any halogen
		},
	},

	"Benzyl_alcohol_deprotection" : 
	{
		"heavy_atoms_diff":7, 
		"rings_diff":1,
		"atomtype":
		{
			"removed" : { "C":7, "H":6 },
		},
	},

}

def check_reaction_types(types):
	for reaction_type in types:
		if reaction_type not in SUPPORTED_CHEMISTRY:
			logger.warning(f"Can't check chemistry of unsupported {reaction_type=}")

def check_chemistry(reaction_type, reactants, product, debug=True):

	assert reaction_type in SUPPORTED_CHEMISTRY
	assert reactants
	assert product	

	CHEMISTRY = SUPPORTED_CHEMISTRY[reaction_type]
	
	if 'heavy_atoms_diff' in CHEMISTRY:
		check = check_count_diff("heavy_atoms", reaction_type, reactants, product, debug=debug)
		if not check:
			return False

	if 'rings_diff' in CHEMISTRY:
		check = check_count_diff("rings", reaction_type, reactants, product, debug=debug)
		if not check:
			return False

	if 'atomtype' in CHEMISTRY:
		check = check_atomtype_diff(reaction_type, reactants, product, debug=debug)
		if not check:
			return False

	if debug:
		logger.success(f'{reaction_type}: All OK')
		
	return True

def check_count_diff(check_type, reaction_type, reactants, product, debug=False):

	# get target value
	diff = SUPPORTED_CHEMISTRY[reaction_type][f'{check_type}_diff']

	# get attribute name
	attr = f'num_{check_type}'

	# get values
	reac_count = getattr(reactants, attr)
	prod_count = getattr(product, attr)
	if debug: logger.var(f'#{check_type} reactants', reac_count)
	if debug: logger.var(f'#{check_type} product', prod_count)

	# check against target value
	if reac_count - diff != prod_count:
		if debug: logger.error(f'{reaction_type}: #{check_type} FAIL')
		return False

	else:
		if debug: logger.success(f'{reaction_type}: #{check_type} OK')

	return True

def check_atomtype_diff(reaction_type, reactants, product, debug=False):

	check_type = 'atomtype'

	# get values
	reac = reactants.atomtype_dict
	prod = product.atomtype_dict

	if debug: 
		logger.var('reactants.atomtype_dict', str(reac))
		logger.var('product.atomtype_dict', str(prod))
	
	if 'removed' in SUPPORTED_CHEMISTRY[reaction_type]['atomtype']:
		removal = check_specific_atomtype_diff(reaction_type, prod, reac, removal=True, debug=debug)

		if not removal:
			return False

	if 'added' in SUPPORTED_CHEMISTRY[reaction_type]['atomtype']:
		addition = check_specific_atomtype_diff(reaction_type, prod, reac, removal=False, debug=debug)

		if not addition:
			return False

	if debug: 
		logger.success(f'{reaction_type}: atomtypes OK')
	
	return True

def check_specific_atomtype_diff(reaction_type, prod, reac, removal=False, debug=False):
	
	if removal:
		add_str = 'removed'
	else:
		add_str = 'added'

	add_dict = SUPPORTED_CHEMISTRY[reaction_type]['atomtype'][add_str]

	if not add_dict:
		return True

	if debug: 
		logger.var(add_str, str(add_dict))
	
	for symbol, count in add_dict.items():

		if symbol == 'Ha':
			p_count = halogen_count(prod)
			r_count = halogen_count(reac)

		else:
			p_count = prod[symbol] if symbol in prod else 0
			r_count = reac[symbol] if symbol in reac else 0
		
		if removal and r_count - p_count != count:
			if debug: 
				logger.error(f'{symbol}: {r_count=} - {p_count=} = {r_count - p_count}')
				logger.error(f'{reaction_type}: atomtype removal {symbol} x {count} FAIL')
			return False

		elif not removal and p_count - r_count != count:
			if debug: 
				logger.error(f'{symbol}: {p_count=} - {r_count=} = {p_count - r_count}')
				logger.error(f'{reaction_type}: atomtype addition {symbol} x {count} FAIL')
			return False
	
	return True

def halogen_count(atomtype_dict):
	count = 0
	symbols = ['F', 'Cl', 'Br', 'I']
	for symbol in symbols:
		if symbol in atomtype_dict:
			count += atomtype_dict[symbol]
	return count

class InvalidChemistryError(Exception):
	pass
