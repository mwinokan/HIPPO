
import mcol

class Reaction:

	def __init__(self,
		db,
		id: int,
		type: str,
		product: int,
		product_amount: float,
	):
		
		self.db = db
		self.id = id
		self.type = type
		self.product = product
		self.product_amount = product_amount
		
	### FACTORIES

	### PROPERTIES

	### METHODS

	### DUNDERS

	def __repr__(self):
		# return f'Reaction(#{self.id})'
		return f'{mcol.bold}{mcol.underline}R{self.id}{mcol.unbold}{mcol.ununderline}'
