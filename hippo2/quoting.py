
import pycule
import mout
import mcol
import json
from pprint import pprint
from pathlib import Path
from requests import TooManyRedirects

from .bb import BuildingBlock, PricePicker

class Quoter:

	def __init__(self, supplier, username, password, force=False):

		self.supplier = supplier
		self.username = username
		self.password = password
		# self.catalogue = catalogue
		self.force = force

		mout.var('supplier', supplier)
		if force:
			mout.warning('Not using cache!')
		mout.var('pycule', str(Path(pycule.__file__).parent), valCol=mcol.file)

		self.requests_data = dict(id_queries={}, smiles_queries={})

		if self.supplier == 'enamine':
			self.query = self.get_bb_info
			self.wrapper = pycule.core.EnamineWrapper(username=username, password=password)

			from pycule.decorators import ENAMINE_MAXIMUM_REQUESTS_PER_MINUTE, ENAMINE_MAXIMUM_REQUESTS_PER_DAY
			mout.var('max frequency', ENAMINE_MAXIMUM_REQUESTS_PER_MINUTE, unit='per minute')
			mout.var('max requests', ENAMINE_MAXIMUM_REQUESTS_PER_DAY, unit='per day')
		
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

	def in_cache(self, comp):
		if comp.name_is_smiles and comp.smiles in self.requests_data['smiles_queries']:
			return True

		elif not comp.name_is_smiles and comp.name in self.requests_data['smiles_queries']:
			return True

		elif not comp.name_is_smiles and comp.name in self.requests_data['id_queries']:
			return True

		return False

	def get_bb_info(self, comp):

		try:

			# ### CACHE

			if comp.name_is_smiles and comp.smiles in self.requests_data['smiles_queries']:
				entry = self.requests_data['smiles_queries'][comp.smiles]
				if entry is None:
					raise NotInCatalogues(f'Cached entry is None: {comp.smiles}')
				comp.name = entry

			elif not comp.name_is_smiles and comp.name in self.requests_data['smiles_queries']:
				entry = self.requests_data['smiles_queries'][comp.name]
				if entry is None:
					raise NotInCatalogues(f'Cached entry is None: {comp.name}')
				comp.name = entry
            
			### QUERY

			# we have a name
			if not comp.name_is_smiles:

				comp_id = comp.name

				if not self.force and comp_id in self.requests_data['id_queries']:
					result = self.requests_data['id_queries'][comp_id]
					# mout.out(f"{mcol.success}using cache.")
					return self.parse_enamine_response(result, comp)

				else:

					# mout.header(f'Quoting: {comp.name}, {comp.smiles}')

					# mout.out(f"Enamine ID search...", end='')
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
						mout.header(f'{mcol.error}not found ({comp_id=}).')
				
			# search by smiles
				
			# mout.out("Enamine SMILES search...", end=' ')

			smiles = comp.smiles
			
			# if not self.force and smiles in self.requests_data['smiles_queries']:
				# result = self.requests_data['smiles_queries'][smiles]

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
				comp_id = self.pick_enamine_exact_data(smiles, result['response']['data'])
				mout.header(f'{mcol.success}found in BB ({smiles} ==> {comp_id}).')
				return comp_id
			except NoDataInReponse:
				pass
		
		result = self.wrapper.exactsearch(smiles, "SCR")

		if result['response']['result']['code'] == 0:
			try:
				comp_id = self.pick_enamine_exact_data(smiles, result['response']['data'])
				mout.header(f'{mcol.success}found in SCR ({smiles} ==> {comp_id}).')
				return comp_id
			except NoDataInReponse:
				pass

		result = self.wrapper.exactsearch(smiles, "MADE")

		if result['response']['result']['code'] == 0:
			try:
				comp_id = self.pick_enamine_exact_data(smiles, result['response']['data'])
				mout.header(f'{mcol.success}found in MADE ({smiles} ==> {comp_id}).')
				return comp_id
			except NoDataInReponse:
				pass

		result = self.wrapper.exactsearch(smiles, "REAL")

		if result['response']['result']['code'] == 0:
			try:
				comp_id = self.pick_enamine_exact_data(smiles, result['response']['data'])
				mout.header(f'{mcol.success}found in REAL ({smiles} ==> {comp_id}).')
				return comp_id
			except NoDataInReponse:
				pass

		mout.header(f'{mcol.error}not found ({smiles=}).')
		return None

	def enamine_comp_id_query(self, compound):

		# try BB
		result = self.wrapper.compoundidsearch(compound.name, catalogue="BB")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.header(f'{mcol.success}found in BB ({compound.name}).')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try SCR
		result = self.wrapper.compoundidsearch(compound.name, catalogue="SCR")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.header(f'{mcol.success}found in SCR ({compound.name}).')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try MADE
		result = self.wrapper.compoundidsearch(compound.name, catalogue="MADE")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.header(f'{mcol.success}found in MADE ({compound.name}).')
				return self.parse_enamine_response(result, compound), result
			except NoDataInReponse:
				pass

		# try REAL
		result = self.wrapper.compoundidsearch(compound.name, catalogue="REAL")
		# pprint(result)

		if result['response']['result']['code'] == 0:
			try:
				mout.header(f'{mcol.success}found in REAL ({compound.name}).')
				return self.parse_enamine_response(result, compound), result
				# return self.parse_enamine_REAL_response(result, compound), result
			except NoDataInReponse:
				pass

		raise NotInCatalogues(f"Couldn't find {compound.name=} in BB, SCR, MADE, or REAL")

	# def parse_enamine_REAL_response(self, result, compound):

	# 	pprint(result)

	# 	delivery = data['deliveryDays']
	# 	days = self.parse_enamine_delivery_string(delivery)
	# 	packs = [dict()]
		
	# 	raise Exception('needs work')
	
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

		if not data['packs'] and 'synthesis' in delivery:
			packs = [dict(amount=1, price=207)]

		else:
			packs = []
			for pack in data['packs']:
				packs.append(self.parse_enamine_ice_pack(pack))

			packs = [p for p in packs if p['price'] > 0]

			if not packs:
				raise NoDataInReponse

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

		if string == '3 weeks synthesis time':
			return 15

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

		# new_data = [d for d in data if not any([p['price'] == 0.0 for p in d['packs']])]

		# if len(new_data) == 0:
		# 	print(data)

		# data = new_data
        
		if len(data) == 0:
			raise NoDataInReponse

		elif len(data) == 1:
			return data[0]

		data = [d for d in data if not any([p['price'] == 0.0 for p in d['packs']])]

		if len(data) > 1:
			mout.warning(f'Taking first data entry after stripping non-priced ({comp_id})')
		return data[0]
		
		# raise Exception(f'ambiguous enamine response {data}')
		
		# assert len(data) == 1, data
		return data[0]

	def pick_enamine_exact_data(self, smiles, data):
	
		if len(data) == 1:
			return data[0]['Id']

		if len(data) == 0:
			raise NoDataInResponse(smiles)

		# sort the results by increasing lead_time, and increasing stereo-complexity
		data = sorted(data, key=lambda x: (self.parse_enamine_delivery_string(x['deliveryDays']), stereo_complexity(x['smile'])))
		mout.warning(f'Picking quickest and least stereo complex result from: {[x["Id"] for x in data]}')
		return data[0]['Id']

def stereo_complexity(smiles):
	return smiles.count('@') + smiles.count('/') + smiles.count('\\')

class NoDataInReponse(Exception):
	...

class NotInCatalogues(Exception):
	...
