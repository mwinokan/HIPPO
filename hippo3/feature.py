
import mcol

from dataclasses import dataclass

from .target import Target

@dataclass
class Feature:

	id: int
	family: Target
	target: int
	chain_name: str
	residue_name: str
	residue_number: int
	atom_names: str

	def __str__(self):
		return f'F{self.id}'

	def __repr__(self):
		return f'{mcol.bold}{mcol.underline}{self.family} {self.chain_name} {self.residue_name} {self.residue_number} [{self.atom_names}]{mcol.unbold}{mcol.ununderline}'

	@property
	def chain_res_name_number_str(self):
		return f'{self.chain_name} {self.residue_name} {self.residue_number}'
