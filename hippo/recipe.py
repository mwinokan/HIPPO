
from dataclasses import dataclass, field

from .compound import Ingredient

import logging
logger = logging.getLogger('HIPPO')

# @dataclass
class Recipe:
	""" A Recipe stores data corresponding to a specific synthetic recipe involving several products, reactants, intermediates, and reactions."""

	def __init__(self, products, reactants, intermediates, reactions):

		from .cset import IngredientSet
		from .rset import ReactionSet
		
		# check typing
		assert isinstance(products, IngredientSet)
		assert isinstance(reactants, IngredientSet)
		assert isinstance(intermediates, IngredientSet)
		assert isinstance(reactions, ReactionSet)
		
		self._products = products
		self._reactants = reactants
		self._intermediates = intermediates
		self._reactions = reactions

	### PROPERTIES

	@property
	def products(self):
		return self._products

	@property
	def reactants(self):
		return self._reactants

	@property
	def intermediates(self):
		return self._intermediates

	@property
	def reactions(self):
		return self._reactions

	@property
	def quotes(self):
		return self.get_quotes()

	@property
	def price(self):
		total = 0
		quotes = self.quotes
		assert len((currencies := set([q.currency for q in quotes]))) == 1, 'Multiple currencies'
		return sum([q.price for q in quotes]), list(currencies)[0]

	@property
	def lead_time(self):
		total = 0
		quotes = self.quotes
		return max([q.lead_time for q in quotes])
	
	### METHODS

	def draw(self, color_mapper = None, node_size=300, graph_only=False):

		"""draw graph of the reaction network"""

		import networkx as nx

		color_mapper = color_mapper or {} 
		colors = {}
		sizes = {}

		graph = nx.DiGraph()

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				key = str(reactant)
				# ingredient = self.get_ingredient(id=reactant.id)
				# key = ingredient.compound_price_amount_str
				graph.add_node(key)

				if not graph_only:
					sizes[key] = self.get_ingredient(id=reactant.id).amount
					if key in color_mapper:
						colors[key] = color_mapper[key]
					else:
						colors[key] = (0.7,0.7,0.7)

		for product in self.products:
			key = str(product)
			# ingredient = self.get_ingredient(id=reactant.id)
			# key = ingredient.compound_price_amount_str
			graph.add_node(key)
			if not graph_only:
				sizes[key] = product.amount
				if key in color_mapper:
					colors[key] = color_mapper[key]
				else:
					colors[key] = (0.7,0.7,0.7)

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				graph.add_edge(str(reactant), str(reaction.product))

		# rescale sizes
		if not graph_only:
			s_min = min(sizes.values())
			sizes = [s/s_min*node_size for s in sizes.values()]

		if graph_only:
			return graph
		else:
			import matplotlib as plt
			# return nx.draw(graph, pos, with_labels=True, font_weight='bold')
			# pos = nx.spring_layout(graph, iterations=200, k=30)
			pos = nx.spring_layout(graph)
			return nx.draw(graph, pos=pos, with_labels=True, font_weight='bold', node_color=list(colors.values()), node_size=sizes)

	def sankey(self):
		"""draw a plotly Sankey diagram"""

		graph = self.draw(graph_only=True)

		import plotly.graph_objects as go

		nodes = {}

		for edge in graph.edges:

			c = edge[0]
			if c not in nodes:
				nodes[c] = len(nodes)

			c = edge[1]
			if c not in nodes:
				nodes[c] = len(nodes)

		source = [nodes[a] for a,b in graph.edges]
		target = [nodes[b] for a,b in graph.edges]
		value = [1 for l in graph.edges]
		
		labels = list(nodes.keys())

		fig = go.Figure(data=[go.Sankey(
				node = dict(
				# pad = 15,
				# thickness = 20,
				# line = dict(color = "black", width = 0.5),
				label = labels,
				# color = "blue"
			),
				link = dict(
				source = source,
				target = target,
				value = value,
		))])

		return fig

	def summary(self):
		"""Print a summary of this recipe"""

		import mcol
		
		logger.header('Recipe')

		logger.var('\n#products', len(self.products))
		for product in self.products:
			logger.var(str(product), f'{product.amount:.2f}', dict(unit='mg'))

		logger.var('\n#intermediates', len(self.intermediates))
		for intermediate in self.intermediates:
			logger.var(str(intermediate), f'{intermediate.amount:.2f}', dict(unit='mg'))

		logger.var('\n#reactants', len(self.reactants))
		for reactant in self.reactants:
			logger.var(str(reactant), f'{reactant.amount:.2f}', dict(unit='mg'))
			# logger.out(f'{mcol.varName}{reactant}{mcol.clear} = {reactant.amount:.2f} {mcol.varType}mg{mcol.clear}, {self.get_reactant_reactions(reactant)}')

		logger.var('\n#reactions', len(self.reactions))
		for reaction in self.reactions:
			logger.var(str(reaction), reaction.reaction_str, dict(unit=reaction.type))

	def get_ingredient(self, id):
		"""Get an ingredient by its compound ID"""
		return [r for r in self.reactants + self.products + self.intermediates if r.id == id][0]

	def get_quotes(self, df=False):

		quotes = [i.quote for i in self.reactants]

		if df:
			from pandas import DataFrame
			return DataFrame([q.dict for q in quotes])

		return quotes

	def write_json(self):

		"""Serialise this recipe object and write it to disk"""

		raise NotImplementedError

	# def get_reactant_reactions(self, reactant):
	# 	"""Get reactions that a reactant is involved in"""
	# 	return [r for r in self.reactions if reactant in r.reactants]

	# def get_reactant_reaction_string(self, reactant):
	# 	reactions = self.get_reactant_reactions(reactant)
	# 	return ' '.join([str(r) for r in reactions])

	### DUNDERS

	def __repr__(self):
		return f'Recipe({self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions})'

		