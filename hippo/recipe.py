
from dataclasses import dataclass, field

from .compound import Ingredient

import logging
logger = logging.getLogger('HIPPO')

# @dataclass
class Recipe:
	""" A Recipe stores data corresponding to a specific synthetic recipe involving several products, reactants, intermediates, and reactions."""

	_db = None

	def __init__(self, db, products, reactants, intermediates, reactions):

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
		self._db = db

	### PROPERTIES

	@property
	def db(self):
		return self._db

	@property
	def products(self):
		return self._products

	@products.setter
	def products(self, a):
		self._products = a

	@property
	def reactants(self):
		return self._reactants

	@reactants.setter
	def reactants(self, a):
		self._reactants = a

	@property
	def intermediates(self):
		return self._intermediates

	@intermediates.setter
	def intermediates(self, a):
		self._intermediates = a

	@property
	def reactions(self):
		return self._reactions

	@reactions.setter
	def reactions(self, a):
		self._reactions = a

	@property
	def quotes(self):
		return self.get_quotes()

	@property
	def price(self):
		total = 0
		quotes = self.quotes
		if not quotes:
			return None
		assert len((currencies := set([q.currency for q in quotes]))) == 1, 'Multiple currencies'
		return sum([q.price for q in quotes]), list(currencies)[0]

	@property
	def lead_time(self):
		total = 0
		quotes = self.quotes
		if not quotes:
			return None
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

				# print(reactant.get_dict(mol=False, metadata=False, poses=False))

				ingredient = self.get_ingredient(id=reactant.id)

				print(ingredient)

				graph.add_node(key, id=reactant.id, smiles=reactant.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time)

				if not graph_only:
					sizes[key] = self.get_ingredient(id=reactant.id).amount
					if key in color_mapper:
						colors[key] = color_mapper[key]
					else:
						colors[key] = (0.7,0.7,0.7)

		for product in self.products:
			key = str(product.compound)
			# ingredient = self.get_ingredient(id=reactant.id)
			# key = ingredient.compound_price_amount_str
			# graph.add_node(key)
			
			ingredient = self.get_ingredient(id=product.id)
			
			print(ingredient)

			graph.add_node(key, id=product.id, smiles=product.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time)

			if not graph_only:
				sizes[key] = product.amount
				if key in color_mapper:
					colors[key] = color_mapper[key]
				else:
					colors[key] = (0.7,0.7,0.7)

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				graph.add_edge(str(reactant), str(reaction.product), reaction_id=reaction.id)

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

	def sankey(self, title=None):
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
		edgedata = [graph.edges[a,b]["reaction_id"] for a,b in graph.edges]

		# print(graph.nodes)
		
		labels = list(nodes.keys())
		
		# compound_ids = [n.id for n in nodes]
		# smiles = [n.smiles for n in nodes]
		# customdata = [(n.id, n.smiles) for n in ]

		hoverkeys = None

		customdata = []
		for key in nodes.keys():
			n = graph.nodes[key]
			
			if not hoverkeys:
				hoverkeys = list(n.keys())
			
			if not n:
				logger.error(f'problem w/ node {key=}')
				compound_id = int(key[1:])
				customdata.append((compound_id, None))

			else:
				# customdata.append((n['id'], n['smiles']))
				d = tuple(v if v is not None else 'N/A' for v in n.values())
				customdata.append(d)
				# id=product.id, smiles=product.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time

		# print(customdata)

		hoverlines = []
		for i,key in enumerate(hoverkeys):
			hoverlines.append(f'{key}=%''{'f'customdata[{i}]''}')
		hovertemplate = 'Compound '+'<br>'.join(hoverlines)+'<extra></extra>'

		# print(hovertemplate)

		# compound_ids = [int(s[1:]) for s in labels]

		# from .cset import CompoundSet
		# smiles = CompoundSet(self.db, compound_ids).smiles

		# print(compound_ids)

		fig = go.Figure(data=[go.Sankey(
				node = dict(
				# pad = 15,
				# thickness = 20,
				# line = dict(color = "black", width = 0.5),
				label = labels,
				# color = "blue"
				customdata=customdata,
				# customdata = ["Long name A1", "Long name A2", "Long name B1", "Long name B2",
                    # "Long name C1", "Long name C2"],
				# hovertemplate='Compound %{label}<br><br>smiles=%{customdata}<extra></extra>',
				hovertemplate=hovertemplate,
			),
				link = dict(
				customdata=edgedata,
				hovertemplate='reaction_id=%{customdata}',
				source = source,
				target = target,
				value = value,
		))])

		if not title:
			title = f"Recipe<br><sup>price={self.price}, lead-time={self.lead_time}</sup>"

		fig.update_layout(title=title)

# link = dict(
#       source = [0, 1, 0, 2, 3, 3], # indices correspond to labels, eg A1, A2, A2, B1, ...
#       target = [2, 3, 3, 4, 4, 5],
#       value = [8, 4, 2, 8, 4, 2],
#       customdata = ["q","r","s","t","u","v"],
#       hovertemplate='Link from node %{source.customdata}<br />'+
#         'to node%{target.customdata}<br />has value %{value}'+
#         '<br />and data %{customdata}<extra></extra>',
#   )

		return fig

	def summary(self):
		"""Print a summary of this recipe"""

		import mcol
		
		logger.header('Recipe')

		price = self.price
		if price:
			logger.var('\nprice', price[0], dict(unit=price[1]))
			logger.var('lead-time', self.lead_time, dict(unit='working days'))

		logger.var('\n#products', len(self.products))
		for product in self.products:
			logger.var(str(product.compound), f'{product.amount:.2f}', dict(unit='mg'))

		logger.var('\n#intermediates', len(self.intermediates))
		for intermediate in self.intermediates:
			logger.var(str(intermediate.compound), f'{intermediate.amount:.2f}', dict(unit='mg'))

		logger.var('\n#reactants', len(self.reactants))
		for reactant in self.reactants:
			logger.var(str(reactant.compound), f'{reactant.amount:.2f}', dict(unit='mg'))
			# logger.out(f'{mcol.varName}{reactant}{mcol.clear} = {reactant.amount:.2f} {mcol.varType}mg{mcol.clear}, {self.get_reactant_reactions(reactant)}')

		logger.var('\n#reactions', len(self.reactions))
		for reaction in self.reactions:
			logger.var(str(reaction), reaction.reaction_str, dict(unit=reaction.type))

	def get_ingredient(self, id):
		"""Get an ingredient by its compound ID"""
		matches = [r for r in self.reactants + self.products + self.intermediates if r.id == id]
		assert len(matches) == 1
		return matches[0]

	def get_quotes(self, df=False):

		quotes = [quote for i in self.reactants if (quote := i.quote)]

		if df:
			from pandas import DataFrame
			return DataFrame([q.dict for q in quotes])

		return quotes

	def write_json(self, file, indent=4):

		"""Serialise this recipe object and write it to disk

		Store
		=====

		- Path to database
		- Timestamp
		- Reactants (& their quotes, amounts)
		- Intermediates (& their quotes)
		- Products (& their poses/scores/fingerprints)
		- Reactions
		- Total Price
		- Lead time

		"""

		import json
		from datetime import datetime

		data = {}

		# Database
		data['database'] = str(self.db.path.resolve())
		data['timestamp'] = str(datetime.now())

		# Recipe properties
		data['price'] = self.price
		data['lead_time'] = self.price

		# IngredientSet's
		data['reactants'] = self.reactants.df.to_json(indent=indent)
		data['intermediates'] = self.intermediates.df.to_json(indent=indent)
		data['products'] = self.products.df.to_json(indent=indent)
		
		# ReactionSet
		data['reaction_ids'] = self.reactions.ids

		logger.writing(file)
		json.dump(open(file, 'wt'), data, indent=indent)

	def copy(self):
		return Recipe(self.db, self.products.copy(), self.reactants.copy(), self.intermediates.copy(), self.reactions.copy())

	# def get_reactant_reactions(self, reactant):
	# 	"""Get reactions that a reactant is involved in"""
	# 	return [r for r in self.reactions if reactant in r.reactants]

	# def get_reactant_reaction_string(self, reactant):
	# 	reactions = self.get_reactant_reactions(reactant)
	# 	return ' '.join([str(r) for r in reactions])

	### DUNDERS

	def __repr__(self):
		return f'Recipe({self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions})'
		# return f'Recipe()'

	def __add__(self, other):
		result = self.copy()
		result.reactants += other.reactants
		result.intermediates += other.intermediates
		result.reactions += other.reactions
		result.products += other.products
		return result

