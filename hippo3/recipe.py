
from dataclasses import dataclass, field

from .compound import Ingredient

import logging
logger = logging.getLogger('HIPPO')

@dataclass
class Recipe:
	""" A Recipe stores data corresponding to a specific synthetic recipe involving one or more reactions."""

	# list of Ingredient options
	products: list = field(default_factory=list)

	# list of Ingredient objects
	reactants: list = field(default_factory=list)
	
	# list of Reaction objects
	reactions: list = field(default_factory=list)

	def draw(self, color_mapper = None):

		"""draw graph of the reaction network"""

		import networkx as nx

		color_mapper = color_mapper or {} 
		colors = {}

		graph = nx.DiGraph()

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				key = str(reactant)
				graph.add_node(key)
				if key in color_mapper:
					colors[key] = color_mapper[key]
				else:
					colors[key] = (0.7,0.7,0.7)

		for product in self.products:
			key = str(product)
			graph.add_node(key)
			if key in color_mapper:
				colors[key] = color_mapper[key]
			else:
				colors[key] = (0.7,0.7,0.7)

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				graph.add_edge(str(reactant), str(reaction.product))

		logger.var('#nodes', len(graph))
		logger.var('#colors', len(colors))

		pos = nx.spring_layout(graph)

		import matplotlib as plt
		# return nx.draw(graph, pos, with_labels=True, font_weight='bold')
		return nx.draw(graph, pos, with_labels=True, font_weight='bold', node_color=list(colors.values()))

	def summary(self):
		
		logger.header('Recipe')

		logger.var('\n#products', len(self.products))
		for product in self.products:
			logger.var(str(product), product.amount, dict(unit='mg'))

		logger.var('\n#reactants', len(self.reactants))
		for reactant in self.reactants:
			logger.var(str(reactant), reactant.amount, dict(unit='mg'))

		logger.var('\n#reactions', len(self.reactions))
		for reaction in self.reactions:
			logger.var(str(reaction), reaction.reaction_str, dict(unit=reaction.type))