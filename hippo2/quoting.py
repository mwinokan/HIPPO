
import pycule
import mout
import json

from .bb import BuildingBlock

class Quoter:

	def __init__(self, supplier, username, password, catalogue, force=False):

		self.supplier = supplier
		self.username = username
		self.password = password
		self.catalogue = catalogue
		self.force = force

		self.requests_data = dict(id_queries={}, smiles_queries={})

		if self.supplier == 'enamine':
			self.query = self.get_bb_info
			self.wrapper = pycule.core.EnamineWrapper(username=username, password=password)
		
		else:
			mout.error(f'Unsupported supplier: "{supplier}"')

	@classmethod
	def from_json(cls, path):

		self = cls.__new__(cls)

		data = json.load(open(path,'rt'))

		self.__init__(
			data['supplier'],
			data['username'],
			data['password'],
			data['catalogue'],
			data['force'],
		)

		self.requests_data = dict(id_queries=data['id_queries'], smiles_queries=data['smiles_queries'])
		
		return self

	def write_json(self, path):
		f = open(path, 'wt')

		data = dict(
			supplier=self.supplier,
			username=self.username,
			password=self.password,
			catalogue=self.catalogue,
			force=self.force,
			id_queries=self.requests_data['id_queries'],
			smiles_queries=self.requests_data['smiles_queries'],
		)

		json.dump(data, f)

	def get_bb_info(self, comp):

		print(comp)

		assert isinstance(comp, BuildingBlock)

		if comp.name_is_smiles:

			print("quote by smiles")

			smiles = comp.smiles

			if not self.force and smiles in self.requests_data['smiles_queries']:
				return self.requests_data['smiles_queries'][smiles]

			result = self.wrapper.similaritysearch(smiles, catalogue=self.catalogue, similarity_value=1.0)
			self.requests_data['smiles_queries'][smiles] = result



			return result

		else:

			print("quote by catalog ID")

			comp_id = comp.name

			if not self.force and comp_id in self.requests_data['id_queries']:
				result = self.requests_data['id_queries'][comp_id]
			else:
				result = self.wrapper.compoundidsearch(comp_id, catalogue=self.catalogue)

			status = parse_enamine_response(result, comp)

			if status:
				self.requests_data['id_queries'][comp_id] = result

			else:
				mout.error(f'Failed to query by ID {comp_id} in "{self.catalogue}"')
				return None
			
			# if result['response']['result']['code']:
			# 	mout.error(f'Failed to query by ID {comp_id} in "{self.catalogue}"')
			# 	return None

			# # get the compound metadata

			# delivery = result['response']['data'][0]['deliveryDays']
			# mout.var("delivery", delivery)

			# delivery = result['response']['data'][0]['deliveryDays']
			# mout.var("delivery", delivery)



			return result

	def __call__(self, *args, **kwargs):
		return self.query(*args, **kwargs)

def parse_enamine_response(result, compound):

	if result['response']['result']['code']:
		mout.error(f'Failed to query by ID {comp_id} in "{self.catalogue}"')
		mout.error(result['response']['result']['message'])
		return False

	delivery = result['response']['data'][0]['deliveryDays']
	mout.var("delivery", delivery, unit='days')

	days = parse_enamine_delivery_string(delivery)
	mout.var("days", days, unit='days')

	packs = []
	for pack in result['response']['data']['packs']:
		packs.append(parse_enamine_ice_pack(pack))

	purity = float(result['response']['data'][0]['purity'])
	mout.var("purity", purity, unit='%')

	return True

def parse_enamine_delivery_string(string):

	if string.startswith('regular delivery, ') and string.endswith('days'):
		return int(string.removeprefix('regular delivery, ').removesuffix(' days'))

	raise Exception

def parse_enamine_ice_pack(pack):

	
	