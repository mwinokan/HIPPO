
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

	# list of Ingredient objects
	intermediates: list = field(default_factory=list)
	
	# list of Reaction objects
	reactions: list = field(default_factory=list)

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
				graph.add_node(key)
				sizes[key] = self.get_ingredient(id=reactant.id).amount
				if key in color_mapper:
					colors[key] = color_mapper[key]
				else:
					colors[key] = (0.7,0.7,0.7)

		for product in self.products:
			key = str(product)
			graph.add_node(key)
			sizes[key] = product.amount
			if key in color_mapper:
				colors[key] = color_mapper[key]
			else:
				colors[key] = (0.7,0.7,0.7)

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				graph.add_edge(str(reactant), str(reaction.product))

		# rescale sizes
		s_min = min(sizes.values())
		sizes = [s/s_min*node_size for s in sizes.values()]

		# pos = nx.spring_layout(graph, iterations=200, k=30)
		pos = nx.spring_layout(graph)

		if graph_only:
			return graph
		else:
			import matplotlib as plt
			# return nx.draw(graph, pos, with_labels=True, font_weight='bold')
			return nx.draw(graph, pos=pos, with_labels=True, font_weight='bold', node_color=list(colors.values()), node_size=sizes)

	def sankey(self):

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
		
		logger.header('Recipe')

		logger.var('\n#products', len(self.products))
		for product in self.products:
			logger.var(str(product), product.amount, dict(unit='mg'))

		logger.var('\n#intermediates', len(self.intermediates))
		for intermediate in self.intermediates:
			logger.var(str(intermediate), intermediate.amount, dict(unit='mg'))

		logger.var('\n#reactants', len(self.reactants))
		for reactant in self.reactants:
			logger.var(str(reactant), reactant.amount, dict(unit='mg'))

		logger.var('\n#reactions', len(self.reactions))
		for reaction in self.reactions:
			logger.var(str(reaction), reaction.reaction_str, dict(unit=reaction.type))

	def get_ingredient(self, id):
		return [r for r in self.reactants + self.products + self.intermediates if r.id == id][0]
