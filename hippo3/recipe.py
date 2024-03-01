
from dataclasses import dataclass, field

from .compound import Ingredient

@dataclass
class Recipe:
	""" A Recipe stores data corresponding to a specific synthetic recipe involving one or more reactions."""

	# the product is a Compound with a specified amount 
	product: Ingredient

	# list of Ingredient objects
	reactants: list = field(default_factory=list)
	
	# list of Reaction objects
	reactions: list = field(default_factory=list)
