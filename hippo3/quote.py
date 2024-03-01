
import mcol
from dataclasses import dataclass, field, asdict

import logging
logger = logging.getLogger('HIPPO')

@dataclass
class Quote:

	compound: int
	smiles: str
	supplier: str
	catalogue: str
	entry: str
	amount: float
	price: float
	currency: str
	purity: float
	lead_time: int
	date: str = field(default_factory=None)

	def __post_init__(self):
		from datetime import datetime
		quote_age = (datetime.today() - datetime.strptime(self.date, '%Y-%m-%d')).days
		if quote_age > 30:
			logger.warning(f'Quote is older than {quote_age} days')
			logger.warning(self)

	def __repr__(self):
		if self.catalogue:
			return f'{mcol.bold}{mcol.underline}{self.supplier}:{self.catalogue}:{self.entry} {self.amount:>8}mg @ {self.purity:.0%} = {self.price:>8} {self.currency} ({self.lead_time} days) {self.smiles}{mcol.unbold}{mcol.ununderline}'
		else:
			return f'{mcol.bold}{mcol.underline}{self.supplier}:{self.entry} {self.amount:>8}mg @ {self.purity:.0%} = {self.price:>8} {self.currency} ({self.lead_time} days) {self.smiles}{mcol.unbold}{mcol.ununderline}'
		
	def asdict(self):
		return asdict(self)
