
from dataclasses import dataclass, field

@dataclass
class Recipe:
	""" A recipe is a set of ingredients and reactions"""

	ingredients: list = field(default_factory=list)
	reactions: list = field(default_factory=list)