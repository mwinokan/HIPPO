
import pycule
from .quote import Quote
import requests

ENAMINE_CATALOGUES = [
	'BB',
	'SCR',
	'MADE',
	'REAL',
]

import logging
logger = logging.getLogger('HIPPO')

class Quoter:

	def __init__(self,
		supplier: str, 
		username: str = None, 
		password: str = None,
		token: str = None,
	):

		if supplier == 'Enamine':

			self.supplier = 'Enamine'
			self.catalogues = ENAMINE_CATALOGUES
			self.get_quote = self.get_enamine_quote

			assert username
			assert password

			self.wrapper = pycule.core.EnamineWrapper(username=username, password=password)

		elif supplier == 'MCule':

			self.supplier = 'MCule'
			self.get_quote = self.get_mcule_quote
			
			assert token

			self.wrapper = pycule.core.MCuleWrapper(authorisation_token=token)
		
		else:
			logger.error(f'Unsupported supplier: "{supplier}"')
			raise NotImplementedError

	### ENAMINE

	def get_enamine_quote(self, compound):

		try:
		
			smiles = compound.smiles

			logger.header(f'Exact search: {smiles=}')

			for catalogue in self.catalogues:

				logger.header(f'Searching in {self.supplier} {catalogue}...')

				try:
					result = self.wrapper.exactsearch(compound.smiles, catalogue)
				except requests.ConnectionError as e:
					logger.error(f'ConnectionError: {e}')
					return None

				if result['response']['result']['code'] != 0:
					continue

				data = result['response']['data']

				assert len(data) == 1

				data = data[0]

				name = data['Id']

				logger.success(f'Found in {catalogue} w/ id={name}')

				break

			else:

				logger.error('No match found')
				return None
			
			### get the information

			for catalogue in self.catalogues:
				
				logger.header(f'Searching by ID in {self.supplier} {catalogue}...')

				result = self.wrapper.compoundidsearch(name, catalogue=catalogue)

				if result['response']['result']['code'] == 0:
					try:
						logger.success(f'Found in {catalogue}')
						return self.parse_enamine_response(result, compound, catalogue), result
					except NoDataInReponse:
						pass

			else:

				logger.error('No catalog entry')
				return None

		except KeyboardInterrupt:
			logger.warning('Interrupted quoting')

	def pick_enamine_exact_data(self, smiles, data):
	
		if len(data) == 1:
			return data[0]['Id']

		if len(data) == 0:
			raise NoDataInResponse(smiles)

		# sort the results by increasing lead_time, and increasing stereo-complexity
		data = sorted(data, key=lambda x: (self.parse_enamine_delivery_string(x['deliveryDays']), stereo_complexity(x['smile'])))
		mout.warning(f'Picking quickest and least stereo complex result from: {[x["Id"] for x in data]}')
		return data[0]['Id']

	def parse_enamine_response(self, result, compound, catalogue):

		quotes = []

		for i,data in enumerate(result['response']['data']):

			
			if i == 0:

				# update compound metadata
				metadata = dict()
				metadata['alias'] = data['name']
				metadata['cas'] = data['cas']
				metadata['mfcd'] = data['mfcd']
				metadata['formula'] = data['formula']
				metadata['molecular_weight'] = data['mw']
				compound.metadata.update(metadata)

				"""
				do something with: 
				* availability
				* storageCond
				* productUrl
				* lastUpdate
				"""
			
			lead_time = self.parse_enamine_delivery_string(data['deliveryDays'])
			entry = data['Id']
			smiles = data['smile']
			purity = data['purity'] / 100
			
			for pack in data['packs']:
				self.parse_enamine_pack(compound, entry, purity, catalogue, pack, lead_time, smiles)

	def pick_enamine_data(self, comp_id,data):
		
		data = [d for d in data if d['Id'] == comp_id]
		
		if len(data) == 0:
			raise NoDataInReponse

		elif len(data) == 1:
			return data[0]

		data = [d for d in data if not any([p['price'] == 0.0 for p in d['packs']])]

		if len(data) > 1:
			logger.warning(f'Taking first data entry after stripping non-priced ({comp_id})')
		return data[0]

	def parse_enamine_delivery_string(self, string):

		if string.startswith('regular delivery, ') and string.endswith('days'):
			return int(string.removeprefix('regular delivery, ').removesuffix(' days'))

		if string.startswith('backorder, ') and string.endswith('days'):
			return int(string.removeprefix('backorder, ').removesuffix(' days'))

		if string == '3 weeks synthesis time':
			return 15

		raise Exception(f'Unexpected delivery string: {string}')

	def parse_enamine_pack(self,compound, entry, purity, catalogue, pack, lead_time, smiles):

		supplier = 'Enamine'

		match pack['measure']:
			case 'g':
				amount = pack['amount'] * 1000
			case 'mg':
				amount = pack['amount']
			case _:
				raise NotImplementedError(f'Not supported pack measure: {pack["measure"]}')

		price = pack['price']
		currency = pack['currencyName']

		if pack['status'] not in ['Normal', 'Ice pack']:
			print(pack)
			raise Exception(f"{pack['status']=}")

		if not price:
			logger.warning(f'Skipping price-less Enamine pack ({entry})')
			return None

		compound.db.insert_quote(
			compound=compound,
			supplier=supplier,
			catalogue=catalogue,
			entry=entry,
			amount=amount,
			price=price,
			currency=currency,
			purity=purity,
			lead_time=lead_time,
			smiles=smiles,
		)

	### MCULE

	def get_mcule_quote(self, compound):

		try:

			# single exact query

			smiles = compound.smiles

			logger.header(f'Exact search: {smiles=}')

			result = self.wrapper.singlequerysearch(smiles)

			if result['response']['results']:

				results = result['response']['results']

				logger.header(results)

				mcule_ids = [data['mcule_id'] for data in results]
				smiles = [data['smiles'] for data in results]

				logger.success(f'Found w/ ids={mcule_ids}')

				for mcule_id,smile in zip(mcule_ids, smiles):

					# logger.header(f'{self.wrapper.compoundavailability(mcule_id)=}')
					# logger.header(f'{self.wrapper.compounddetails(mcule_id)=}')
					result = self.wrapper.compoundprices(mcule_id)

					for pack in result['response']['best_prices']:
						self.parse_mcule_pack(compound, mcule_id, pack, smile)

					# logger.header(f'{self.wrapper.compoundpricesamount(mcule_id)=}')

			return result

		except KeyboardInterrupt:
			logger.warning('Interrupted quoting')

	def parse_mcule_pack(self, compound, entry, pack, smiles):

		supplier = 'MCule'

		match pack['unit']:
			case 'g':
				amount = pack['amount'] * 1000
			case 'mg':
				amount = pack['amount']
			case _:
				raise NotImplementedError(f'Not supported pack measure: {pack["unit"]}')

		price = pack['price']
		currency = pack['currency']
		purity = pack['purity'] / 100
		lead_time = pack['delivery_time_working_days']

		compound.db.insert_quote(
			compound=compound,
			supplier=supplier,
			catalogue=None,
			entry=entry,
			amount=amount,
			price=price,
			currency=currency,
			purity=purity,
			lead_time=lead_time,
			smiles=smiles,
		)

class NoDataInReponse(Exception):
	...

class NotInCatalogues(Exception):
	...
