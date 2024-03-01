
import mcol
from dataclasses import dataclass

@dataclass
class Target:

	id: int
	name: str

	def __str__(self):
		return f'T{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self} "{self.name}"{mcol.unbold}{mcol.ununderline}'
