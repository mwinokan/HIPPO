from dataclasses import dataclass, field

from .compound import Ingredient

import mcol

from tqdm import tqdm

import logging

logger = logging.getLogger("HIPPO")


# @dataclass
class Recipe:
    """A Recipe stores data corresponding to a specific synthetic recipe involving several products, reactants, intermediates, and reactions."""

    _db = None

    def __init__(self, db, *, products, reactants, intermediates, reactions):

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
        self._hash = None

        self._score = None

        # caches
        self._product_compounds = None
        self._product_poses = None
        self._product_interactions = None

    ### FACTORIES

    @classmethod
    def from_reaction(
        cls,
        reaction,
        amount=1,
        *,
        debug: bool = False,
        pick_cheapest: bool = True,
        permitted_reactions=None,
        quoted_only: bool = False,
        supplier: None | str = None,
        unavailable_reaction="error",
        reaction_checking_cache=None,
        reaction_reactant_cache=None,
        inner=False,
    ):
        """

        :param reaction:
        :param amount:  (Default value = 1)
        :param *:
        :param debug: bool:  (Default value = False)
        :param pick_cheapest: bool:  (Default value = True)
        :param permitted_reactions:  (Default value = None)
        :param quoted_only: bool:  (Default value = False)
        :param supplier: None | str:  (Default value = None)
        :param unavailable_reaction:  (Default value = 'error')
        :param reaction_checking_cache:  (Default value = None)
        :param reaction_reactant_cache:  (Default value = None)
        :param inner:  (Default value = False)

        """

        from .reaction import Reaction

        assert isinstance(reaction, Reaction)

        from .cset import IngredientSet
        from .rset import ReactionSet

        if debug:
            logger.debug(f"Recipe.from_reaction({amount=}, {pick_cheapest=})")
        if debug:
            logger.debug(f"{reaction.id=}")
        if debug:
            logger.debug(f"{reaction.product.id=}")
        if debug:
            logger.debug(f"{amount=}")
        if debug:
            logger.debug(f"{reaction.reactants.ids=}")

        if permitted_reactions:
            assert reaction in permitted_reactions
            # raise NotImplementedError

        db = reaction.db

        recipe = cls.__new__(cls)
        recipe.__init__(
            db,
            products=IngredientSet(db, [reaction.product.as_ingredient(amount=amount)]),
            reactants=IngredientSet(db, [], supplier=supplier),
            intermediates=IngredientSet(db, []),
            reactions=ReactionSet(db, [reaction.id], sort=False),
        )

        recipes = [recipe]

        if quoted_only or supplier:
            if debug:
                logger.debug(f"Checking reactant_availability: {reaction=}")
            if reaction_checking_cache and reaction.id in reaction_checking_cache:
                ok = reaction_checking_cache[reaction.id]
                print("reaction_checking_cache used")
            else:
                ok = reaction.check_reactant_availability(supplier=supplier)
                # print('cache not used')
                if reaction_checking_cache is not None:
                    reaction_checking_cache[reaction.id] = ok
            if not ok:
                if unavailable_reaction == "error":
                    logger.error(f"Reactants not available for {reaction=}")
                if pick_cheapest:
                    return None
                else:
                    return []

        def get_reactant_amount_pairs(reaction):
            """

            :param reaction:

            """
            if reaction_reactant_cache and reaction.id in reaction_reactant_cache:
                print("reaction_reactant_cache used")
                return reaction_reactant_cache[reaction.id]
            else:
                pairs = reaction.get_reactant_amount_pairs(compound_object=False)
                if reaction_reactant_cache is not None:
                    reaction_reactant_cache[reaction.id] = pairs
                return pairs

        # logger.debug(f'get_reactant_amount_pairs({reaction.id})')
        # pairs = reaction.get_reactant_amount_pairs()
        pairs = get_reactant_amount_pairs(reaction)

        for reactant, reactant_amount in pairs:

            reactant = db.get_compound(id=reactant)

            if debug:
                logger.debug(f"{reactant.id=}, {reactant_amount=}")

            # scale amount
            reactant_amount *= amount
            reactant_amount /= reaction.product_yield

            inner_reactions = reactant.get_reactions(
                none="quiet", permitted_reactions=permitted_reactions
            )

            if inner_reactions:

                if debug:
                    if len(inner_reactions) == 1:
                        logger.debug(f"Reactant has ONE inner reaction")
                    else:
                        logger.warning(f"{reactant=} has MULTIPLE inner reactions")

                new_recipes = []

                inner_recipes = []
                for reaction in inner_reactions:
                    reaction_recipes = Recipe.from_reaction(
                        reaction=reaction,
                        amount=reactant_amount,
                        debug=debug,
                        pick_cheapest=False,
                        quoted_only=quoted_only,
                        supplier=supplier,
                        unavailable_reaction=unavailable_reaction,
                        reaction_checking_cache=reaction_checking_cache,
                        reaction_reactant_cache=reaction_reactant_cache,
                        inner=True,
                    )
                    inner_recipes += reaction_recipes

                for recipe in recipes:

                    for inner_recipe in inner_recipes:

                        combined_recipe = recipe.copy()

                        combined_recipe.reactants += inner_recipe.reactants
                        combined_recipe.intermediates += inner_recipe.intermediates
                        combined_recipe.reactions += inner_recipe.reactions
                        combined_recipe.intermediates.add(
                            reactant.as_ingredient(reactant_amount, supplier=supplier)
                        )

                        new_recipes.append(combined_recipe)

                recipes = new_recipes

            else:

                ingredient = reactant.as_ingredient(reactant_amount, supplier=supplier)
                for recipe in recipes:
                    recipe.reactants.add(ingredient)

        # reverse ReactionSet's
        if not inner:
            for recipe in recipes:
                recipe.reactions.reverse()

        if pick_cheapest:
            priced = [r for r in recipes if r.get_price(supplier=supplier)]
            # priced = [r for r in recipes if r.price]
            if not priced:
                logger.error("0 recipes with prices, can't choose cheapest")
                return recipes
            return sorted(priced, key=lambda r: r.get_price(supplier=supplier))[0]
            # return sorted(priced, key=lambda r: r.price)[0]

        return recipes

    @classmethod
    def from_reactions(
        cls,
        reactions,
        amount=1,
        pick_cheapest: bool = True,
        permitted_reactions=None,
        final_products_only=True,
        return_products=False,
        debug=False,
    ):
        """

        :param reactions:
        :param amount:  (Default value = 1)
        :param pick_cheapest: bool:  (Default value = True)
        :param permitted_reactions:  (Default value = None)
        :param final_products_only:  (Default value = True)
        :param return_products:  (Default value = False)
        :param debug:  (Default value = False)

        """

        from .rset import ReactionSet
        from .cset import IngredientSet, CompoundSet

        assert isinstance(reactions, ReactionSet)

        db = reactions.db

        if debug:
            logger.debug("Recipe.from_reactions()")
            logger.var("reactions", reactions)
            logger.var("amount", amount)
            logger.var("final_products_only", final_products_only)
            logger.var("permitted_reactions", permitted_reactions)

        # get all the products
        products = reactions.products

        if debug:
            logger.var("products", products)

        # return products

        if final_products_only:

            if debug:
                logger.var("products.str_ids", products.str_ids)

            # raise NotImplementedError
            ids = reactions.db.execute(
                f"""
                SELECT DISTINCT compound_id FROM compound
                LEFT JOIN reactant ON compound_id = reactant_compound
                WHERE reactant_compound IS NULL
                AND compound_id IN {products.str_ids}
            """
            ).fetchall()

            ids = [i for i, in ids]

            products = CompoundSet(db, ids)
            if debug:
                logger.var("final products", products)

            # return ids

            if return_products:
                return products

        recipe = Recipe.from_compounds(
            compounds=products,
            amount=amount,
            permitted_reactions=reactions,
            pick_cheapest=pick_cheapest,
        )

        return recipe

    @classmethod
    def from_compounds(
        cls,
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
        unavailable_reaction="error",
        reaction_checking_cache=None,
        reaction_reactant_cache=None,
    ):
        """

        :param compounds:
        :param amount: float:  (Default value = 1)
        :param debug: bool:  (Default value = False)
        :param pick_cheapest: bool:  (Default value = True)
        :param permitted_reactions:  (Default value = None)
        :param quoted_only: bool:  (Default value = False)
        :param supplier: None | str:  (Default value = None)
        :param solve_combinations: bool:  (Default value = True)
        :param pick_first: bool:  (Default value = False)
        :param warn_multiple_solutions: bool:  (Default value = True)
        :param pick_cheapest_inner_routes: bool:  (Default value = False)
        :param unavailable_reaction:  (Default value = 'error')
        :param reaction_checking_cache:  (Default value = None)
        :param reaction_reactant_cache:  (Default value = None)

        """

        from tqdm import tqdm

        from .cset import CompoundSet

        assert isinstance(compounds, CompoundSet)

        # if permitted_reactions:
        #   raise NotImplementedError

        n_comps = len(compounds)

        assert n_comps

        if not hasattr(amount, "__iter__"):
            amount = [amount] * n_comps

        options = []

        if debug:
            logger.var("#compounds", n_comps)
            logger.info("Solving individual compound recipes...")

        if n_comps > 1:
            generator = tqdm(zip(compounds, amount), total=n_comps)
        else:
            generator = zip(compounds, amount)

        for comp, a in generator:

            comp_options = []

            for reaction in comp.reactions:

                if permitted_reactions and reaction not in permitted_reactions:
                    continue

                sol = Recipe.from_reaction(
                    reaction=reaction,
                    amount=a,
                    pick_cheapest=pick_cheapest_inner_routes,
                    debug=debug,
                    permitted_reactions=permitted_reactions,
                    quoted_only=quoted_only,
                    supplier=supplier,
                    unavailable_reaction=unavailable_reaction,
                    reaction_checking_cache=reaction_checking_cache,
                    reaction_reactant_cache=reaction_reactant_cache,
                )

                if pick_cheapest_inner_routes:
                    if sol:
                        comp_options.append(sol)
                else:
                    assert isinstance(sol, list)
                    comp_options += sol

            if warn_multiple_solutions and len(comp_options) > 1:
                logger.warning(f"Multiple solutions for compound={comp}")

            if not comp_options:
                logger.error(f"No solutions for compound={comp}")
            else:
                options.append(comp_options)

        assert all(options)

        from itertools import product

        if debug:
            logger.info("Solving recipe combinations...")
        combinations = list(product(*options))

        if not solve_combinations:
            return combinations

        if pick_first:
            combinations = [combinations[0]]

        if debug:
            logger.info("Combining recipes...")

        solutions = []

        if n_comps > 1:
            generator = tqdm(combinations)
        else:
            generator = combinations

        for combo in generator:

            if debug:
                logger.debug(f"Combination of {len(combo)} recipes")

            solution = combo[0]

            for recipe in combo[1:]:
                solution += recipe

            solutions.append(solution)

        if not solutions:
            logger.error("No solutions!")

        if pick_first:
            return solutions[0]

        if pick_cheapest:
            if debug:
                logger.info("Picking cheapest...")
            priced = [r for r in solutions if r.price]
            if not priced:
                logger.error("0 recipes with prices, can't choose cheapest")
                return solutions
            return sorted(priced, key=lambda r: r.price)[0]

        return solutions

    @classmethod
    def from_reactants(cls, reactants, amount=1, debug=False, return_products=False):
        """

        :param reactants:
        :param amount:  (Default value = 1)
        :param debug:  (Default value = False)
        :param return_products:  (Default value = False)

        """

        from .cset import IngredientSet

        if isinstance(reactants, IngredientSet):
            reactant_ids = reactants.compound_ids
        else:
            reactant_ids = reactants.ids

        db = reactants.db

        all_reactants = set(reactant_ids)

        possible_reactions = []

        # recursively search for possible reactions
        for i in range(300):

            if debug:
                logger.debug(i)

            # reaction_ids = db.get_possible_reaction_ids(compound_ids=compound_ids)
            reaction_ids = db.get_possible_reaction_ids(compound_ids=all_reactants)

            if not reaction_ids:
                break

            if debug:
                logger.debug(f"Adding {len(reaction_ids)} reactions")

            possible_reactions += reaction_ids

            if debug:
                logger.var("reaction_ids", reaction_ids)

            product_ids = db.get_possible_reaction_product_ids(
                reaction_ids=reaction_ids
            )

            if debug:
                logger.var("product_ids", product_ids)

            n_prev = len(all_reactants)

            all_reactants |= set(product_ids)

            if n_prev == len(all_reactants):
                break

        else:
            raise NotImplementedError("Maximum recursion depth exceeded")

        possible_reactions = list(set(possible_reactions))

        if debug:
            logger.var("all possible reactions", possible_reactions)

        from .rset import ReactionSet

        rset = ReactionSet(db, possible_reactions, sort=False)

        recipe = cls.from_reactions(
            rset,
            amount=amount,
            permitted_reactions=rset,
            pick_cheapest=False,
            debug=debug,
            return_products=return_products,
        )

        return recipe

    @classmethod
    def from_ingredients(cls, db, ingredients):
        """

        :param db:
        :param ingredients:

        """
        raise NotImplementedError
        self = cls.__new__(cls)
        self.__init__(...)
        return self

    @classmethod
    def from_json(
        cls,
        db,
        path,
        debug=True,
        allow_db_mismatch=False,
        clear_quotes=False,
        data=None,
        db_mismatch_warning: bool = True,
    ):
        """

        :param db:
        :param path:
        :param debug:  (Default value = True)
        :param allow_db_mismatch:  (Default value = False)
        :param clear_quotes:  (Default value = False)
        :param data:  (Default value = None)

        """

        # imports
        import json
        from .cset import IngredientSet
        from .rset import ReactionSet

        # load JSON
        if not data:
            if debug:
                logger.reading(path)
            data = json.load(open(path, "rt"))

        # check metadata
        if str(db.path.resolve()) != data["database"]:
            if db_mismatch_warning:
                logger.var("session", str(db.path.resolve()))
                logger.var("in file", data["database"])
            if allow_db_mismatch:
                if db_mismatch_warning:
                    logger.warning("Database path mismatch")
            else:
                logger.error(
                    "Database path mismatch, set allow_db_mismatch=True to ignore"
                )
                return None

        if debug:
            logger.info(f'Recipe was generated at: {data["timestamp"]}')
        price = data["price"]

        # IngredientSets
        products = IngredientSet.from_ingredient_dicts(db, data["products"])
        intermediates = IngredientSet.from_ingredient_dicts(db, data["intermediates"])
        reactants = IngredientSet.from_ingredient_dicts(
            db, data["reactants"], supplier=data["reactant_supplier"]
        )

        if clear_quotes:
            reactants.df["quote_id"] = None
            reactants.df["quoted_amount"] = None

        # ReactionSet
        reactions = ReactionSet(db, data["reaction_ids"], sort=False)

        if debug:
            logger.var("reactants", reactants)
            logger.var("intermediates", intermediates)
            logger.var("products", products)
            logger.var("reactions", reactions)

        # Create the object
        self = cls.__new__(cls)
        self.__init__(
            db,
            products=products,
            reactants=reactants,
            intermediates=intermediates,
            reactions=reactions,
        )

        return self

    ### PROPERTIES

    @property
    def db(self):
        """ """
        return self._db

    @property
    def products(self):
        """ """
        return self._products

    @property
    def product_poses(self) -> "PoseSet":
        if self._product_poses is None:
            self._product_poses = self.product_compounds.poses
        return self._product_poses

    @property
    def product_compounds(self) -> "CompoundSet":
        if self._product_compounds is None:
            self._product_compounds = self.products.compounds
        return self._product_compounds

    @property
    def product_interactions(self) -> "InteractionSet":
        """Product pose interactions"""
        if self._product_interactions is None:
            self._product_interactions = self.product_poses.interactions
        return self._product_interactions

    @property
    def product(self):
        """ """
        assert len(self.products) == 1
        return self.products[0]

    @products.setter
    def products(self, a):
        """

        :param a:

        """
        self._products = a
        self.__flag_modification()

    @property
    def reactants(self):
        """ """
        return self._reactants

    @reactants.setter
    def reactants(self, a):
        """

        :param a:

        """
        self._reactants = a
        self.__flag_modification()

    @property
    def intermediates(self):
        """ """
        return self._intermediates

    @intermediates.setter
    def intermediates(self, a):
        """

        :param a:

        """
        self._intermediates = a
        self.__flag_modification()

    @property
    def reactions(self):
        """ """
        return self._reactions

    @reactions.setter
    def reactions(self, a):
        """

        :param a:

        """
        self._reactions = a
        self.__flag_modification()

    @property
    def price(self):
        """ """
        # total = 0
        # quotes = self.quotes
        # if not quotes:
        #   return None
        # assert len((currencies := set([q.currency for q in quotes]))) == 1, 'Multiple currencies'
        # return sum([q.price for q in quotes]), list(currencies)[0]
        return self.reactants.price

    @property
    def num_products(self):
        return len(self.products)

    @property
    def hash(self):
        return self._hash

    @property
    def score(self):
        return self._score

    # @property
    # def lead_time(self):
    #   total = 0
    #   quotes = self.quotes
    #   if not quotes:
    #       return None
    #   return max([q.lead_time for q in quotes])

    ### METHODS

    def get_price(self, supplier=None):
        """

        :param supplier:  (Default value = None)

        """
        return self.reactants.get_price(supplier=supplier)

    def draw(self, color_mapper=None, node_size=300, graph_only=False):
        """draw graph of the reaction network

        :param color_mapper:  (Default value = None)
        :param node_size:  (Default value = 300)
        :param graph_only:  (Default value = False)

        """

        import networkx as nx

        color_mapper = color_mapper or {}
        colors = {}
        sizes = {}

        graph = nx.DiGraph()

        for reaction in self.reactions:
            for reactant in reaction.reactants:
                key = str(reactant)
                ingredient = self.get_ingredient(id=reactant.id)

                graph.add_node(
                    key,
                    id=reactant.id,
                    smiles=reactant.smiles,
                    amount=ingredient.amount,
                    price=ingredient.price,
                    lead_time=ingredient.lead_time,
                )

                if not graph_only:
                    sizes[key] = self.get_ingredient(id=reactant.id).amount
                    if key in color_mapper:
                        colors[key] = color_mapper[key]
                    else:
                        colors[key] = (0.7, 0.7, 0.7)

        for product in self.products:
            key = str(product.compound)
            ingredient = self.get_ingredient(id=product.id)

            graph.add_node(
                key,
                id=product.id,
                smiles=product.smiles,
                amount=ingredient.amount,
                price=ingredient.price,
                lead_time=ingredient.lead_time,
            )

            if not graph_only:
                sizes[key] = product.amount
                if key in color_mapper:
                    colors[key] = color_mapper[key]
                else:
                    colors[key] = (0.7, 0.7, 0.7)

        for reaction in self.reactions:
            for reactant in reaction.reactants:
                graph.add_edge(
                    str(reactant),
                    str(reaction.product),
                    id=reaction.id,
                    type=reaction.type,
                    product_yield=reaction.product_yield,
                )

        # rescale sizes
        if not graph_only:
            s_min = min(sizes.values())
            sizes = [s / s_min * node_size for s in sizes.values()]

        if graph_only:
            return graph
        else:
            import matplotlib as plt

            # return nx.draw(graph, pos, with_labels=True, font_weight='bold')
            # pos = nx.spring_layout(graph, iterations=200, k=30)
            pos = nx.spring_layout(graph)
            return nx.draw(
                graph,
                pos=pos,
                with_labels=True,
                font_weight="bold",
                node_color=list(colors.values()),
                node_size=sizes,
            )

    def sankey(self, title=None):
        """draw a plotly Sankey diagram

        :param title:  (Default value = None)

        """

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

        source = [nodes[a] for a, b in graph.edges]
        target = [nodes[b] for a, b in graph.edges]
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
                logger.error(f"problem w/ node {key=}")
                compound_id = int(key[1:])
                customdata.append((compound_id, None))

            else:
                # customdata.append((n['id'], n['smiles']))
                d = tuple(v if v is not None else "N/A" for v in n.values())
                customdata.append(d)
                # id=product.id, smiles=product.smiles, amount=ingredient.amount, price=ingredient.price, lead_time=ingredient.lead_time

        hoverkeys_edges = None

        # edgedata = [graph.edges[a,b]["reaction_id"] for a,b in graph.edges]

        customdata_edges = []

        for s, t in graph.edges.keys():
            edge = graph.edges[s, t]

            if not hoverkeys_edges:
                hoverkeys_edges = list(edge.keys())

            if not n:
                logger.error(f"problem w/ edge {s=} {t=}")
                customdata_edges.append((None, None, None))

            else:
                d = tuple(v if v is not None else "N/A" for v in edge.values())
                customdata_edges.append(d)

        hoverlines = []
        for i, key in enumerate(hoverkeys):
            hoverlines.append(f"{key}=%" "{" f"customdata[{i}]" "}")
        hovertemplate = "Compound " + "<br>".join(hoverlines) + "<extra></extra>"

        hoverlines_edges = []
        for i, key in enumerate(hoverkeys_edges):
            hoverlines_edges.append(f"{key}=%" "{" f"customdata[{i}]" "}")
        hovertemplate_edges = (
            "Reaction " + "<br>".join(hoverlines_edges) + "<extra></extra>"
        )

        # print(hovertemplate)

        # compound_ids = [int(s[1:]) for s in labels]

        # from .cset import CompoundSet
        # smiles = CompoundSet(self.db, compound_ids).smiles

        # print(compound_ids)

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(
                        # pad = 15,
                        # thickness = 20,
                        # line = dict(color = "black", width = 0.5),
                        label=labels,
                        # color = "blue"
                        customdata=customdata,
                        # customdata = ["Long name A1", "Long name A2", "Long name B1", "Long name B2",
                        # "Long name C1", "Long name C2"],
                        # hovertemplate='Compound %{label}<br><br>smiles=%{customdata}<extra></extra>',
                        hovertemplate=hovertemplate,
                    ),
                    link=dict(
                        customdata=customdata_edges,
                        hovertemplate=hovertemplate_edges,
                        source=source,
                        target=target,
                        value=value,
                    ),
                )
            ]
        )

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

    def summary(self, price: bool = True):
        """Print a summary of this recipe

        :param price: bool:  (Default value = True)

        """

        import mcol

        logger.header("Recipe")

        if price:
            price = self.price
            if price:
                logger.var("\nprice", price[0], dict(unit=price[1]))
                # logger.var('lead-time', self.lead_time, dict(unit='working days'))

        logger.var("\n#products", len(self.products))
        for product in self.products:
            logger.var(str(product.compound), f"{product.amount:.2f}", dict(unit="mg"))

        logger.var("\n#intermediates", len(self.intermediates))
        for intermediate in self.intermediates:
            logger.var(
                str(intermediate.compound),
                f"{intermediate.amount:.2f}",
                dict(unit="mg"),
            )

        logger.var("\n#reactants", len(self.reactants))
        for reactant in self.reactants:
            logger.var(
                str(reactant.compound), f"{reactant.amount:.2f}", dict(unit="mg")
            )
            # logger.out(f'{mcol.varName}{reactant}{mcol.clear} = {reactant.amount:.2f} {mcol.varType}mg{mcol.clear}, {self.get_reactant_reactions(reactant)}')

        logger.var("\n#reactions", len(self.reactions))
        for reaction in self.reactions:
            logger.var(str(reaction), reaction.reaction_str, dict(unit=reaction.type))

    def get_ingredient(self, id):
        """Get an ingredient by its compound ID

        :param id:

        """
        matches = [r for r in self.reactants if r.id == id]
        if not matches:
            matches = [r for r in self.intermediates if r.id == id]
        if not matches:
            matches = [r for r in self.products if r.id == id]

        assert len(matches) == 1
        return matches[0]

    def add_to_all_reactants(self, amount=20):
        """

        :param amount:  (Default value = 20)

        """
        self.reactants.df["amount"] += amount

    def write_json(self, file, *, indent="\t", **kwargs):
        r"""Serialise this recipe object and write it to disk

        :param file:
        :param indent:  (Default value = '\t')

        """
        import json

        data = self.get_dict(serialise_price=True, **kwargs)
        logger.writing(file)
        json.dump(data, open(file, "wt"), indent=indent)

    def get_dict(
        self,
        *,
        price: bool = True,
        reactant_supplier: bool = True,
        database: bool = True,
        timestamp: bool = True,
        compound_ids_only: bool = False,
        products: bool = True,
        serialise_price: bool = False,
    ):
        """Serialise this recipe object

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

        :param *:
        :param price: bool:  (Default value = True)
        :param reactant_supplier: bool:  (Default value = True)
        :param database: bool:  (Default value = True)
        :param timestamp: bool:  (Default value = True)
        :param compound_ids_only: bool:  (Default value = False)
        :param products: bool:  (Default value = True)
        :param serialise_price: bool:  (Default value = False)

        """

        import json
        from datetime import datetime

        data = {}

        # Database
        if database:
            data["database"] = str(self.db.path.resolve())
        if timestamp:
            data["timestamp"] = str(datetime.now())

        # Recipe properties
        if price and serialise_price:
            data["price"] = self.price.get_dict()
        elif price:
            data["price"] = self.price

        if reactant_supplier:
            data["reactant_supplier"] = self.reactants.supplier
        # data['lead_time'] = self.lead_time

        # IngredientSets
        if compound_ids_only:
            data["reactant_ids"] = self.reactants.compound_ids
            data["intermediate_ids"] = self.intermediates.compound_ids
            if products:
                data["products_ids"] = self.products.compound_ids

        else:
            data["reactants"] = self.reactants.df.to_dict(orient="list")
            data["intermediates"] = self.intermediates.df.to_dict(orient="list")
            if products:
                data["products"] = self.products.df.to_dict(orient="list")

        # ReactionSet
        data["reaction_ids"] = self.reactions.ids

        return data

    def write_CAR_csv(self, file, return_df=False):
        """List of reactions for CAR

        Columns:

        * target-name
        * no-steps
        * concentration = None
        * amount-required
        * batch-tag

        per reaction

        * reactant-1-1
        * reactant-2-1
        * reaction-product-smiles-1
        * reaction-name-1
        * reaction-recipe-1
        * reaction-groupby-column-1

        :param file:
        :param return_df:  (Default value = False)

        """

        from .cset import CompoundSet
        from pandas import DataFrame
        from tqdm import tqdm

        # solve each product's reaction

        rows = []

        for product in tqdm(self.products):

            prod_cset = CompoundSet(self.db, [product.id])

            sub_recipes = prod_cset.get_recipes(permitted_reactions=self.reactions)

            for sub_recipe in sub_recipes:

                row = {
                    "target-names": str(product.compound),
                    "no-steps": 0,
                    "concentration-required-mM": None,
                    "amount-required-uL": None,
                    "batch-tag": None,
                }

                for i, reaction in enumerate(sub_recipe.reactions):

                    i = i + 1

                    row["no-steps"] += 1

                    match len(reaction.reactants):
                        case 1:
                            row[f"reactant-1-{i}"] = reaction.reactants[0].smiles
                            row[f"reactant-2-{i}"] = None
                        case 2:
                            row[f"reactant-1-{i}"] = reaction.reactants[0].smiles
                            row[f"reactant-2-{i}"] = reaction.reactants[1].smiles
                        case _:
                            raise NotImplementedError(
                                f"Unsupported number of reactants for {reaction=}: {len(reaction.reactants)}"
                            )

                    row[f"reaction-product-smiles-{i}"] = reaction.product.smiles
                    row[f"reaction-name-{i}"] = reaction.type
                    row[f"reaction-recipe-{i}"] = None
                    row[f"reaction-groupby-column-{i}"] = None
                    # row[f'reaction-id-{i}'] = int(reaction.id)

                rows.append(row)

        df = DataFrame(rows)

        df = df.convert_dtypes()

        for n_steps in set(df["no-steps"]):
            subset = df[df["no-steps"] == n_steps]
            this_file = file.replace(".csv", f"_{n_steps}steps.csv")
            logger.writing(this_file)
            subset.to_csv(this_file, index=False)

        logger.writing(file)
        df.to_csv(file, index=False)

        return df

    def copy(self):
        """ """
        return Recipe(
            self.db,
            products=self.products.copy(),
            reactants=self.reactants.copy(),
            intermediates=self.intermediates.copy(),
            reactions=self.reactions.copy(),
            # supplier=self.supplier
        )

    def get_product_fingerprint(self):
        """Calculate the combined fingerprint of all product poses in this set"""

        poses = self.products.compounds.poses

        null_count = self.db.count_where(
            table="pose", key=f"pose_id IN {poses.str_ids} AND pose_fingerprint IS NULL"
        )

        if null_count:
            logger.warning(f"{null_count} poses have no fingerprint")

        return poses.fingerprint

    def __flag_modification(self):
        self._product_interactions = None
        self._score = None
        self._product_compounds = None
        self._product_poses = None

    ### DUNDERS

    def __repr__(self):
        if self.score:
            return f"Recipe_{self.hash}({self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions}, score={self.score:.3f})"
        elif self.hash:
            return f"Recipe_{self.hash}({self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions})"
        else:
            return f"Recipe({self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions})"
        # return f'Recipe()'

    def __add__(self, other):
        result = self.copy()
        result.reactants += other.reactants
        result.intermediates += other.intermediates
        result.reactions += other.reactions
        result.products += other.products
        return result


class Route(Recipe):
    """A recipe with a single product, that is stored in the database"""

    def __init__(self, db, *, route_id, product, reactants, intermediates, reactions):

        from .cset import IngredientSet
        from .rset import ReactionSet

        # check typing
        assert isinstance(product, IngredientSet)
        assert isinstance(reactants, IngredientSet)
        assert isinstance(intermediates, IngredientSet)
        assert isinstance(reactions, ReactionSet)

        assert len(product) == 1
        assert isinstance(route_id, int)
        assert route_id

        self._id = route_id
        self._products = product
        self._product_id = product.ids[0]
        self._reactants = reactants
        self._intermediates = intermediates
        self._reactions = reactions
        self._db = db

    ### FACTORIES

    @classmethod
    def from_json(cls, db, path, data=None):
        """

        :param db:
        :param path:
        :param data:  (Default value = None)

        """

        import json
        from .cset import IngredientSet
        from .rset import ReactionSet

        if data is None:
            data = json.load(open(path, "rt"))

        self = cls.__new__(cls)

        self._db = db
        self._id = data["id"]

        self._product_id = data["product_id"]
        self._products = IngredientSet.from_compounds(
            compounds=None, ids=[self._product_id], db=db
        )  # IngredientSet

        self._reactants = IngredientSet.from_json(
            db=db,
            path=None,
            data=data["reactants"]["data"],
            supplier=data["reactants"]["supplier"],
        )
        self._intermediates = IngredientSet.from_json(
            db=db,
            path=None,
            data=data["intermediates"]["data"],
            supplier=data["intermediates"]["supplier"],
        )
        self._reactions = ReactionSet(
            db=db, indices=data["reactions"]["indices"]
        )  # ReactionSet

        return self

    ### PROPERTIES

    @property
    def product(self):
        """ """
        return self._products[0]

    @property
    def product_compound(self):
        """ """
        return self.product.compound

    @property
    def id(self):
        """ """
        return self._id

    ### METHODS

    def get_dict(self):
        """Serialisable dictionary"""
        data = {}

        data["id"] = self.id
        data["product_id"] = self.product.id
        data["reactants"] = self.reactants.get_dict()
        data["intermediates"] = self.intermediates.get_dict()
        data["reactions"] = self.reactions.get_dict()

        return data

    ### DUNDERS

    def __repr__(self):
        return f"{mcol.bold}{mcol.underline}Route #{self.id}: {repr(self.product_compound)}"  #' {self.reactants} --> {self.intermediates} --> {self.product_compound} via {self.reactions})'


class RouteSet:
    """A set of Route objects"""

    def __init__(self, db, routes):

        data = {}
        for route in routes:
            assert isinstance(route, Route)
            data[route.id] = route

        self._data = data
        self._db = db

    ### FACTORIES

    @classmethod
    def from_json(cls, db, path, data=None):
        """

        :param db:
        :param path:
        :param data:  (Default value = None)

        """

        self = cls.__new__(cls)

        if data is None:
            import json

            data = json.load(open(path, "rt"))

        new_data = {}
        for d in tqdm(data["routes"].values()):
            route_id = d["id"]
            new_data[route_id] = Route.from_json(db=db, path=None, data=d)

        self._data = new_data
        self._db = db

        return self

    ### PROPERTIES

    @property
    def data(self):
        """ """
        return self._data

    @property
    def db(self):
        """ """
        return self._db

    @property
    def routes(self):
        """ """
        return self.data.values()

    @property
    def product_ids(self):
        """Get the :class:`.Compound` ID's of the products"""
        ids = self.db.select_where(
            table="route",
            query="route_product",
            key=f"route_id IN {self.str_ids}",
            multiple=True,
        )
        return [i for i, in ids]

    @property
    def str_ids(self):
        """Return an SQL formatted tuple string of the :class:`.Route` ID's"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def ids(self):
        """Return the :class:`.Route` IDs"""
        return self.data.keys()

    ### METHODS

    def copy(self):
        """ """
        return RouteSet(self.db, self.data.values())

    def set_db_pointers(self, db):
        """

        :param db:

        """
        self._db = db
        for route in self.data.values():
            route._db = db

    def clear_db_pointers(self):
        """ """
        self._db = None
        for route in self.data.values():
            route._db = None

    def get_dict(self):
        """Get serialisable dictionary"""

        data = dict(db=str(self.db), routes={})

        # populate with routes
        for route_id, route in self.data.items():
            data["routes"][route_id] = route.get_dict()

        return data

    def pop_id(self):
        """ """

        route_id, route = self.data.popitem()

        return route_id

    def pop(self):
        """ """

        route_id, route = self.data.popitem()

        return route

    def shuffle(self):
        """Randomly shuffle the routes in this set"""
        import random

        items = list(self.data.items())
        random.shuffle(items)
        self._data = dict(items)

    ### DUNDERS

    def __len__(self):
        return len(self.data)


class RecipeSet:
    """
    RecipeSet class

    :param param1: this is a first param
    :param param2: this is a second param
    :returns: this is a description of what is returned
    :raises keyError: raises an exception
    """

    def __init__(self, db, directory, pattern="*.json"):

        from pathlib import Path

        self._db = db
        self._json_directory = Path(directory)
        self._json_pattern = pattern

        self._json_paths = {}
        for path in self._json_directory.glob(self._json_pattern):
            self._json_paths[
                path.name.removeprefix("Recipe_").removesuffix(".json")
            ] = path.resolve()

        logger.reading(f"{directory}/{pattern}")

        self._recipes = {}
        for key, path in tqdm(self._json_paths.items()):
            recipe = Recipe.from_json(
                db=self.db,
                path=path,
                allow_db_mismatch=True,
                debug=False,
                db_mismatch_warning=False,
            )
            recipe._hash = key
            self._recipes[key] = recipe

        # print(self._recipes)

    ### FACTORIES

    ### PROPERTIES

    @property
    def db(self):
        return self._db

    ### METHODS

    def get_values(
        self,
        key: str,
        progress: bool = False,
        serialise_price: bool = False,
    ):

        values = []
        recipes = self._recipes.values()

        if progress:
            recipes = tqdm(recipes)

        for recipe in recipes:
            value = getattr(recipe, key)
            if serialise_price and key == "price":
                value = value.amount
            values.append(value)

        return values

    def get_df(self, **kwargs) -> "pandas.DataFrame":

        data = []

        for recipe in self:

            d = recipe.get_dict(
                # reactant_supplier=False,
                database=False,
                timestamp=False,
                **kwargs,
                # timestamp=False,
            )

            data.append(d)

        from pandas import DataFrame

        return DataFrame(data)

    def items(self):
        return self._recipes.items()

    ### DUNDERS

    def __len__(self) -> int:
        """Number of recipes in this set"""
        return len(self._recipes)

    def __getitem__(
        self,
        key: int | str,
    ) -> Recipe:

        match key:

            case int():
                return list(self._recipes.values())[key]

            case str():
                return self._recipes[key]

            case _:
                logger.error(
                    f"Unsupported type for RecipeSet.__getitem__(): {key=} {type(key)}"
                )

        return None

    def __iter__(self):
        return iter(self._recipes.values())
