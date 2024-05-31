
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

	### FACTORIES


	@classmethod
	def from_reaction(cls, 
		reaction, 
		amount=1, 
		debug: bool = False, 
		pick_cheapest: bool = True,
		permitted_reactions=None,
		quoted_only: bool = False,
		supplier: None | str = None,
	):

		from .reaction import Reaction
		assert isinstance(reaction, Reaction)

		from .cset import IngredientSet
		from .rset import ReactionSet

		if debug: logger.debug(f'Recipe.from_reaction({amount=}, {pick_cheapest=})')
		if debug: logger.debug(f'{reaction.id=}')
		if debug: logger.debug(f'{reaction.product.id=}')
		if debug: logger.debug(f'{amount=}')
		if debug: logger.debug(f'{reaction.reactants.ids=}')

		if permitted_reactions:
			assert reaction in permitted_reactions
			# raise NotImplementedError

		db = reaction.db

		recipe = cls.__new__(cls)
		recipe.__init__(db, 
			products=IngredientSet(db, [reaction.product.as_ingredient(amount=amount)]),
			reactants=IngredientSet(db, []),
			intermediates=IngredientSet(db, []),
			reactions=ReactionSet(db, [reaction.id]),
		)

		recipes = [recipe]

		if quoted_only:
			ok = reaction.check_reactant_availability(supplier=supplier)
			if not ok:
				logger.error(f'Reactants not available for {reaction=}')
				if pick_cheapest:
					return None
				else:
					return []
		
		pairs = reaction.get_reactant_amount_pairs()
		for reactant, reactant_amount in pairs:
			if debug: logger.debug(f'{reactant.id=}, {reactant_amount=}')

			# scale amount
			reactant_amount *= amount
			reactant_amount /= reaction.product_yield

			inner_reactions = reactant.get_reactions(none='quiet', permitted_reactions=permitted_reactions)

			if inner_reactions:

				if debug: 
					if len(inner_reactions) == 1:
						logger.debug(f'Reactant has ONE inner reaction')
					else:
						logger.warning(f'{reactant=} has MULTIPLE inner reactions')

				new_recipes = []

				inner_recipes = []
				for reaction in inner_reactions:
					reaction_recipes = Recipe.from_reaction(reaction=reaction, amount=reactant_amount, debug=debug, pick_cheapest=False, quoted_only=quoted_only)
					inner_recipes += reaction_recipes

				for recipe in recipes:

					for inner_recipe in inner_recipes:

						combined_recipe = recipe.copy()

						combined_recipe.reactants += inner_recipe.reactants
						combined_recipe.intermediates += inner_recipe.intermediates
						combined_recipe.reactions += inner_recipe.reactions
						combined_recipe.intermediates.add(reactant.as_ingredient(reactant_amount))

						new_recipes.append(combined_recipe)
					
				recipes = new_recipes

			else:

				ingredient = reactant.as_ingredient(reactant_amount)
				for recipe in recipes:
					recipe.reactants.add(ingredient)

		if pick_cheapest:
			priced = [r for r in recipes if r.get_price(supplier=supplier)]
			if not priced:
				logger.error("0 recipes with prices, can't choose cheapest")
				return recipes
			return sorted(priced, key=lambda r: r.get_price(supplier=supplier))[0]

		return recipes

	@classmethod
	def from_reactions(cls, 
		reactions, 
		amount=1,
		pick_cheapest: bool = True,
		permitted_reactions=None,
		final_products_only=True,
	):

		from .rset import ReactionSet
		from .cset import IngredientSet
		assert isinstance(reactions, ReactionSet)

		# get all the products
		products = reactions.products

		return products

		if final_products_only:
			ids = reactions.db.execute(f'''
				SELECT DISTINCT compound_id FROM compound
				LEFT JOIN reactant ON compound_id = reactant_compound
				WHERE reactant_compound IS NULL
				AND compound_id IN {products.str_ids}
			''').fetchall()

			return ids

		recipe = Recipe.from_compounds(compounds=products, amount=amount, permitted_reactions=reactions, pick_cheapest=pick_cheapest)

		return recipe

	@classmethod
	def from_compounds(cls,
		compounds, 
		amount: float = 1, 
		debug: bool = False, 
		pick_cheapest: bool = True,
		permitted_reactions=None,
		quoted_only: bool = False,
		supplier: None | str = None,
		solve_combinations: bool = True,
		pick_first: bool = False,
		warn_multiple_solutions: bool = True,
		pick_cheapest_inner_routes: bool = False,
	):

		from tqdm import tqdm

		from .cset import CompoundSet
		assert isinstance(compounds, CompoundSet)

		# if permitted_reactions:
		# 	raise NotImplementedError

		n_comps = len(compounds)

		if not hasattr(amount, '__iter__'):
			amount = [amount] * n_comps

		options = []

		logger.debug(f'{n_comps=}')

		logger.info('Solving individual compound recipes...')
		for comp, a in tqdm(zip(compounds, amount), total=n_comps):

			comp_options = []

			for reaction in comp.reactions:

				if permitted_reactions and reaction not in permitted_reactions:
					continue

				if pick_cheapest_inner_routes:
					sol = Recipe.from_reaction(reaction=reaction, amount=a, pick_cheapest=True, debug=debug, permitted_reactions=permitted_reactions, quoted_only=quoted_only, supplier=supplier)
					if sol:
					# assert isinstance(sol, Recipe)
						comp_options.append(sol)
				else:
					sols = Recipe.from_reaction(reaction=reaction, amount=a, pick_cheapest=False, debug=debug, permitted_reactions=permitted_reactions, quoted_only=quoted_only, supplier=supplier)	
					assert isinstance(sols, list)
					comp_options += sols

			if warn_multiple_solutions and len(comp_options) > 1:
				logger.warning(f'Multiple solutions for compound={comp}')

			if not comp_options:
				logger.error(f'No solutions for compound={comp}')
			else:
				options.append(comp_options)

		assert all(options)

		from itertools import product
		logger.info('Solving recipe combinations...')
		combinations = list(product(*options))

		if not solve_combinations:
			return combinations

		if pick_first:
			combinations = [combinations[0]]
			
		logger.info('Combining recipes...')
		
		solutions = []
		for combo in tqdm(combinations):

			logger.debug(f'Combination of {len(combo)} recipes')

			solution = combo[0]

			for recipe in combo[1:]:
				solution += recipe

			solutions.append(solution)

		if not solutions:
			logger.error('No solutions!')

		if pick_first:
			return solutions[0]

		if pick_cheapest:
			logger.info('Picking cheapest...')
			priced = [r for r in solutions if r.price]
			if not priced:
				logger.error("0 recipes with prices, can't choose cheapest")
				return solutions
			return sorted(priced, key=lambda r: r.price)[0]

		return solutions

	@classmethod
	def from_reactants(cls, reactants, amount=1, debug=False):

		# print(reactants.ids)


		from .cset import IngredientSet

		if isinstance(reactants, IngredientSet):
			reactant_ids = reactants.compound_ids
		else:
			reactant_ids = reactants.ids

		db = reactants.db

		compound_ids = reactant_ids

		possible_reactions = []

		# recursively search for possible reactions
		for i in range(300):
			
			reaction_ids = db.get_possible_reaction_ids(compound_ids=compound_ids)

			if not reaction_ids:
				break

			possible_reactions += reaction_ids

			if debug:
				logger.var('reaction_ids', reaction_ids)

			product_ids = db.get_possible_reaction_product_ids(reaction_ids=reaction_ids)

			if debug:
				logger.var('product_ids', product_ids)

			compound_ids = product_ids

		else:
			raise NotImplementedError('Maximum recursion depth exceeded')

		if debug:
			logger.var('all possible reactions', possible_reactions)

		from .rset import ReactionSet

		rset = ReactionSet(db, possible_reactions)

		return cls.from_reactions(rset, amount=amount, permitted_reactions=rset, pick_cheapest=False)

	@classmethod
	def from_ingredients(cls, db, ingredients):
		raise NotImplementedError
		self = cls.__new__(cls)
		self.__init__(...)
		return self

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

	# @property
	# def quotes(self):
		# return self.get_quotes()

	@property
	def price(self):
		# total = 0
		# quotes = self.quotes
		# if not quotes:
		# 	return None
		# assert len((currencies := set([q.currency for q in quotes]))) == 1, 'Multiple currencies'
		# return sum([q.price for q in quotes]), list(currencies)[0]
		return self.reactants.price

	# @property
	# def lead_time(self):
	# 	total = 0
	# 	quotes = self.quotes
	# 	if not quotes:
	# 		return None
	# 	return max([q.lead_time for q in quotes])
	
	### METHODS

	def get_price(self, supplier=None):
		return self.reactants.get_price(supplier=supplier)

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
				ingredient = self.get_ingredient(id=reactant.id)

				graph.add_node(key, id=reactant.id, smiles=reactant.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time)

				if not graph_only:
					sizes[key] = self.get_ingredient(id=reactant.id).amount
					if key in color_mapper:
						colors[key] = color_mapper[key]
					else:
						colors[key] = (0.7,0.7,0.7)

		for product in self.products:
			key = str(product.compound)			
			ingredient = self.get_ingredient(id=product.id)

			graph.add_node(key, id=product.id, smiles=product.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time)

			if not graph_only:
				sizes[key] = product.amount
				if key in color_mapper:
					colors[key] = color_mapper[key]
				else:
					colors[key] = (0.7,0.7,0.7)

		for reaction in self.reactions:
			for reactant in reaction.reactants:
				graph.add_edge(str(reactant), str(reaction.product), id=reaction.id, type=reaction.type, product_yield=reaction.product_yield)

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

		hoverkeys_edges = None
		
		# edgedata = [graph.edges[a,b]["reaction_id"] for a,b in graph.edges]

		customdata_edges = []

		for s,t in graph.edges.keys():
			edge = graph.edges[s,t]
			
			if not hoverkeys_edges:
				hoverkeys_edges = list(edge.keys())
			
			if not n:
				logger.error(f'problem w/ edge {s=} {t=}')
				customdata_edges.append((None, None, None))

			else:
				d = tuple(v if v is not None else 'N/A' for v in edge.values())
				customdata_edges.append(d)

		hoverlines = []
		for i,key in enumerate(hoverkeys):
			hoverlines.append(f'{key}=%''{'f'customdata[{i}]''}')
		hovertemplate = 'Compound '+'<br>'.join(hoverlines)+'<extra></extra>'

		hoverlines_edges = []
		for i,key in enumerate(hoverkeys_edges):
			hoverlines_edges.append(f'{key}=%''{'f'customdata[{i}]''}')
		hovertemplate_edges = 'Reaction '+'<br>'.join(hoverlines_edges)+'<extra></extra>'

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
				customdata=customdata_edges,
				hovertemplate=hovertemplate_edges,
				source = source,
				target = target,
				value = value,
		))])

		if not title:
			# title = f"Recipe<br><sup>price={self.price}, lead-time={self.lead_time}</sup>"
			title = f"Recipe<br><sup>price={self.price}</sup>"

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
			# logger.var('lead-time', self.lead_time, dict(unit='working days'))

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

	def add_to_all_reactants(self, amount=20):
		self.reactants.df['amount'] += amount

	# def get_quotes(self, df=False):

	# 	quotes = [quote for i in self.reactants if (quote := i.quote)]

	# 	if df:
	# 		from pandas import DataFrame
	# 		return DataFrame([q.dict for q in quotes])

	# 	return quotes

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
		# data['lead_time'] = self.lead_time

		# IngredientSets
		data['reactants'] = self.reactants.df.to_dict(orient='records')
		data['intermediates'] = self.intermediates.df.to_dict(orient='records')
		data['products'] = self.products.df.to_dict(orient='records')
		
		# ReactionSet
		data['reaction_ids'] = self.reactions.ids

		# return data

		logger.writing(file)
		json.dump(data, open(file, 'wt'), indent=indent)

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

