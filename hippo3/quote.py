
import mcol
from dataclasses import dataclass

@dataclass
class Quote:

	compound: int
	supplier: str
	catalogue: str
	entry: str
	amount: float
	price: float
	currency: str
	lead_time: int
	date: str

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self.amount} mg of "{self.entry}" in {self.supplier} {self.catalogue} for {self.price} {self.currency} in {self.lead_time} days{mcol.unbold}{mcol.ununderline}'