
from dataclasses import dataclass

@dataclass
class Feature:

	id: int
	type: Target
	target: str
	chain_name: str
	residue_name: str
	residue_number: int
	atom_name: str
	atom_number: int
