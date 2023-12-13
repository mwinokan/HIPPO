
import pycule
import mout
import mcol
import json
from pprint import pprint
from requests import TooManyRedirects

from .bb import BuildingBlock, PricePicker

class Quoter:

	def __init__(self, supplier, username, password, force=False):

		self.supplier = supplier
		self.username = username
		self.password = password
		# self.catalogue = catalogue
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
			data['force'],
		)

		self.requests_data = dict(id_queries=data['id_queries'], smiles_queries=data['smiles_queries'])
		
		return self

	def write_json(self, path):

		mout.out(f'writing {mcol.file}{path}{mcol.clear} ... ', end='')

		f = open(path, 'wt')

		data = dict(
			supplier=self.supplier,
			username=self.username,
			password=self.password,
			# catalogue=self.catalogue,
			force=self.force,
			id_queries=self.requests_data['id_queries'],
			smiles_queries=self.requests_data['smiles_queries'],
		)

		json.dump(data, f)
		mout.out('Done.')

	def get_bb_info(self, comp):

		mout.header(f'Quoting: {comp.name}, {comp.smiles}')

		try:

			# assert isinstance(comp, BuildingBlock), type(BuildingBlock)

			if comp.name_is_smiles and comp.smiles in self.requests_data['smiles_queries']:
				entry = self.requests_data['smiles_queries'][comp.smiles]
				if entry is None:
					raise NotInCatalogues('cached')
				comp.name = entry

			elif not comp.name_is_smiles and comp.name in self.requests_data['smiles_queries']:
				entry = self.requests_data['smiles_queries'][comp.name]
				if entry is None:
					raise NotInCatalogues('cached')
				comp.name = entry
			
			# we have a name
			if not comp.name_is_smiles:

				mout.out(f"Enamine ID search...", end='')

				comp_id = comp.name

				if not self.force and comp_id in self.requests_data['id_queries']:
					result = self.requests_data['id_queries'][comp_id]
					mout.out(f"{mcol.success}using cache.")
					return self.parse_enamine_response(result, comp)

				else:

					try:
						parsed, result = self.enamine_comp_id_query(comp)
					
						if parsed:
							self.requests_data['id_queries'][comp_id] = result
						else:
							mout.error(f'Failed to parse request data {comp_id}')
							return None

						# mout.out(f"{mcol.success}OK.")
						return parsed

					except NotInCatalogues as e:
						# mout.error(f'{e}: {comp.name}')
						mout.out(f'{mcol.error}not found.')
				
			# search by smiles
				
			mout.out("Enamine SMILES search...", end=' ')

			smiles = comp.smiles
			
			if not self.force and smiles in self.requests_data['smiles_queries']:
				result = self.requests_data['smiles_queries'][smiles]

			comp_id = self.enamine_exact_search(comp)

			if comp_id is None:
				self.requests_data['smiles_queries'][comp.name] = None
				raise NotInCatalogues("Couldn't find SMILES in BB, SCR, MADE, or REAL")

			if comp.name_is_smiles:
				self.requests_data['smiles_queries'][smiles] = comp_id
			else:
				self.requests_data['smiles_queries'][comp.name] = comp_id
			
			# mout.out(f"{mcol.success}OK.")

			mout.out(f"Renaming and retrying...")
			comp.name = comp_id

			return self.get_bb_info(comp)

		except TooManyRedirects as e:
			mout.error(f'TooManyRedirects {e}')
			mout.error(f'{comp}')
			return None

	def __call__(self, *args, **kwargs):
		return self.query(*args, **kwargs)

	def enamine_exact_search(self, compound):

		'''
		exactsearch(self, smiles: str, currency: str = 'USD', catalogue: str = 'BB') -> requests.models.Response
        This search returns the exact match search results for the queried compound.
        
        Args:
            code (str): Required -  the SMILES to search for
            currency (str):  Required -  set to “USD”. Available values: "USD" and "EUR"
            mode (str): Optional (default is "BB"). Available values: "BB", "SCR", "REAL"
        Returns:
            dict: dictionary containing the search response
		'''

		smiles = compound.smiles

		result = self.wrapper.exactsearch(smiles, "BB")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				assert len(result['response']['data']) == 1
				mout.out(f'{mcol.success}found in BB.')
				return result['response']['data'][0]['Id']
			except NoDataInReponse:
				pass
		
		result = self.wrapper.exactsearch(smiles, "SCR")

		if result['response']['result']['code'] == 0:
			try:
				assert len(result['response']['data']) == 1
				mout.out(f'{mcol.success}found in SCR.')
				return result['response']['data'][0]['Id']
			except NoDataInReponse:
				pass

		result = self.wrapper.exactsearch(smiles, "MADE")

		if result['response']['result']['code'] == 0:
			try:
				assert len(result['response']['data']) == 1
				mout.out(f'{mcol.success}found in MADE.')
				return result['response']['data'][0]['Id']
			except NoDataInReponse:
				pass

		result = self.wrapper.exactsearch(smiles, "REAL")

		if result['response']['result']['code'] == 0:
			try:
				assert len(result['response']['data']) == 1
				mout.out(f'{mcol.success}found in REAL.')
				return result['response']['data'][0]['Id']
			except NoDataInReponse:
				pass

		mout.out(f'{mcol.error}not found.')
		return None

	def enamine_comp_id_query(self, compound):

		# try BB
		result = self.wrapper.compoundidsearch(compound.name, catalogue="BB")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.out(f'{mcol.success}found in BB.')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try SCR
		result = self.wrapper.compoundidsearch(compound.name, catalogue="SCR")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.out(f'{mcol.success}found in SCR.')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try MADE
		result = self.wrapper.compoundidsearch(compound.name, catalogue="MADE")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.out(f'{mcol.success}found in MADE.')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try REAL
		result = self.wrapper.compoundidsearch(compound.name, catalogue="REAL")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.out(f'{mcol.success}found in REAL.')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		raise NotInCatalogues("Couldn't find comp_id in BB, SCR, MADE, or REAL")

	def parse_enamine_response(self, result, compound):

		from pprint import pprint

		# if result['response']['result']['code']:
		# 	mout.error(f'Failed to query by ID {compound.name}')
		# 	mout.error(result['response']['result']['message'])
		# 	return None

		if len(result['response']['data']) == 1:
			data = result['response']['data'][0]
		else:
			data = self.pick_enamine_data(compound.name, result['response']['data'])
		
		delivery = data['deliveryDays']
		days = self.parse_enamine_delivery_string(delivery)
		# mout.var("days", days, unit='days')

		# purity = float(data['purity'])
		# mout.var("purity", purity, unit='%')
		
		packs = []
		for pack in data['packs']:
			packs.append(self.parse_enamine_ice_pack(pack))
		# pprint(packs)

		#### DO STUFF TO THE COMPOUND

		compound.lead_time = days
		compound.price_picker = PricePicker(packs)

		return dict(days=days, packs=packs)
		# return dict(days=days, purity=purity, packs=packs)

	def parse_enamine_delivery_string(self, string):

		if string.startswith('regular delivery, ') and string.endswith('days'):
			return int(string.removeprefix('regular delivery, ').removesuffix(' days'))

		if string.startswith('backorder, ') and string.endswith('days'):
			return int(string.removeprefix('backorder, ').removesuffix(' days'))

		raise Exception(f'Unexpected delivery string: {string}')

	def parse_enamine_ice_pack(self, pack):

		if pack['measure'] == 'g':
			amount = float(pack['amount'])*1000 # mg
		
		elif pack['measure'] == 'mg':
			amount = float(pack['amount']) # mg

		else:
			raise Exception(f'Unsupported pack_measure: {pack["measure"]}')

		price = float(pack['price']) # USD
		assert pack['currencyName'] == 'USD'

		return dict(amount=amount, price=price)
		
	def pick_enamine_data(self, comp_id,data):
		data = [d for d in data if d['Id'] == comp_id]

		if len(data) == 0:
			raise NoDataInReponse

		elif len(data) == 1:
			return data[0]
		
		raise Exception(f'ambiguous enamine response {data}')
		
		# assert len(data) == 1, data
		return data[0]

	# def is_enamine_data_ok(self, response):

	# 	if result['response']['result']['code']:
	# 		return False

class NoDataInReponse(Exception):
	...

class NotInCatalogues(Exception):
	...
