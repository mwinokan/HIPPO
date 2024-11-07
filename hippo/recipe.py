from dataclasses import dataclass, field

from .compound import Ingredient

import mcol

import mrich


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
        get_ingredient_quotes: bool = True,
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
            mrich.debug(
                f"Recipe.from_reaction(R{reaction.id}, {amount=}, {pick_cheapest=})"
            )
            mrich.debug(f"{reaction.product.id=}")
            mrich.debug(f"{reaction.reactants.ids=}")

        if permitted_reactions:
            assert reaction in permitted_reactions
            # raise NotImplementedError

        db = reaction.db

        recipe = cls.__new__(cls)
        recipe.__init__(
            db,
            products=IngredientSet(
                db,
                [
                    reaction.product.as_ingredient(
                        amount=amount, get_quote=get_ingredient_quotes
                    )
                ],
            ),
            reactants=IngredientSet(db, [], supplier=supplier),
            intermediates=IngredientSet(db, []),
            reactions=ReactionSet(db, [reaction.id], sort=False),
        )

        recipes = [recipe]

        if quoted_only or supplier:
            if debug:
                mrich.debug(f"Checking reactant_availability: {reaction=}")
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
                    mrich.error(f"Reactants not available for {reaction=}")
                if pick_cheapest:
                    return None
                else:
                    return []

        def get_reactant_amount_pairs(reaction):
            if reaction_reactant_cache and reaction.id in reaction_reactant_cache:
                print("reaction_reactant_cache used")
                return reaction_reactant_cache[reaction.id]
            else:
                pairs = reaction.get_reactant_amount_pairs(compound_object=False)
                if reaction_reactant_cache is not None:
                    reaction_reactant_cache[reaction.id] = pairs
                return pairs

        if debug:
            mrich.debug(f"get_reactant_amount_pairs({reaction.id})")
        pairs = get_reactant_amount_pairs(reaction)

        for reactant, reactant_amount in pairs:

            reactant = db.get_compound(id=reactant)

            if debug:
                mrich.debug(f"{reactant.id=}, {reactant_amount=}")

            # scale amount
            reactant_amount *= amount
            reactant_amount /= reaction.product_yield

            inner_reactions = reactant.get_reactions(
                none="quiet", permitted_reactions=permitted_reactions
            )

            if inner_reactions:

                if debug:
                    if len(inner_reactions) == 1:
                        mrich.debug(f"Reactant has ONE inner reaction")
                    else:
                        mrich.warning(f"{reactant=} has MULTIPLE inner reactions")

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
            if debug:
                mrich.debug("Picking cheapest")
            priced = [r for r in recipes if r.get_price(supplier=supplier)]
            # priced = [r for r in recipes if r.price]
            if not priced:
                mrich.error("0 recipes with prices, can't choose cheapest")
                return recipes
            sorted_recipes = sorted(
                priced, key=lambda r: r.get_price(supplier=supplier)
            )

            if debug:
                for recipe in recipes:
                    mrich.debug(f"{recipe}, {recipe.price}")

            return sorted_recipes[0]
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
            mrich.debug("Recipe.from_reactions()")
            mrich.var("reactions", reactions)
            mrich.var("amount", amount)
            mrich.var("final_products_only", final_products_only)
            mrich.var("permitted_reactions", permitted_reactions)

        # get all the products
        products = reactions.products

        if debug:
            mrich.var("products", products)

        # return products

        if final_products_only:

            if debug:
                mrich.var("products.str_ids", products.str_ids)

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
                mrich.var("final products", products)

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
        **kwargs,
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

        from .cset import CompoundSet

        assert isinstance(compounds, CompoundSet)

        # if permitted_reactions:
        #   raise NotImplementedError

        n_comps = len(compounds)

        assert n_comps

        if not hasattr(amount, "__iter__"):
            amount = [amount] * n_comps

        options = []

        mrich.var("#compounds", n_comps)

        if n_comps > 1:
            generator = mrich.track(
                zip(compounds, amount), prefix="Solving individual compound recipes..."
            )
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
                    **kwargs,
                )

                if pick_cheapest_inner_routes:
                    if sol:
                        comp_options.append(sol)
                else:
                    assert isinstance(sol, list)
                    comp_options += sol

            if not comp_options:
                mrich.error(f"No solutions for compound={comp} ({comp.reactions.ids=})")
                continue

            if pick_cheapest and len(comp_options) > 1:
                if warn_multiple_solutions:
                    mrich.warning(f"Multiple solutions for compound={comp}")
                if debug:
                    mrich.debug("Picking cheapest...")
                priced = [r for r in comp_options if r.price]
                comp_options = sorted(priced, key=lambda r: r.price)[:1]

            if warn_multiple_solutions and len(comp_options) > 1:
                mrich.warning(f"Multiple solutions for compound={comp}")
                if debug:
                    mrich.debug(f"{comp_options=}")
            else:
                mrich.success(f"Found solution for compound={comp}")

            options.append(comp_options)

        assert all(options)

        from itertools import product

        mrich.print("Solving recipe combinations...")
        combinations = list(product(*options))

        if not solve_combinations:
            return combinations

        if pick_first:
            combinations = [combinations[0]]

        mrich.print()

        solutions = []

        if n_comps > 1:
            generator = mrich.track(combinations, prefix="Combining recipes...")
        else:
            generator = combinations

        for combo in generator:

            if debug:
                mrich.debug(f"Combination of {len(combo)} recipes")

            solution = combo[0]

            for i, recipe in enumerate(combo[1:]):
                if debug:
                    mrich.debug(i + 1)
                solution += recipe

            solutions.append(solution)

        if not solutions:
            mrich.error("No solutions!")

        if pick_first:
            return solutions[0]

        if pick_cheapest:
            mrich.print("Picking cheapest...")
            priced = [r for r in solutions if r.price]
            if not priced:
                mrich.error("0 recipes with prices, can't choose cheapest")
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
                mrich.debug(i)

            # reaction_ids = db.get_possible_reaction_ids(compound_ids=compound_ids)
            reaction_ids = db.get_possible_reaction_ids(compound_ids=all_reactants)

            if not reaction_ids:
                break

            if debug:
                mrich.debug(f"Adding {len(reaction_ids)} reactions")

            possible_reactions += reaction_ids

            if debug:
                mrich.var("reaction_ids", reaction_ids)

            product_ids = db.get_possible_reaction_product_ids(
                reaction_ids=reaction_ids
            )

            if debug:
                mrich.var("product_ids", product_ids)

            n_prev = len(all_reactants)

            all_reactants |= set(product_ids)

            if n_prev == len(all_reactants):
                break

        else:
            raise NotImplementedError("Maximum recursion depth exceeded")

        possible_reactions = list(set(possible_reactions))

        if debug:
            mrich.var("all possible reactions", possible_reactions)

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
                mrich.reading(path)
            data = json.load(open(path, "rt"))

        # check metadata
        if str(db.path.resolve()) != data["database"]:
            if db_mismatch_warning:
                mrich.var("session", str(db.path.resolve()))
                mrich.var("in file", data["database"])
            if allow_db_mismatch:
                if db_mismatch_warning:
                    mrich.warning("Database path mismatch")
            else:
                mrich.error(
                    "Database path mismatch, set allow_db_mismatch=True to ignore"
                )
                return None

        if debug:
            mrich.print(f'Recipe was generated at: {data["timestamp"]}')
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
            mrich.var("reactants", reactants)
            mrich.var("intermediates", intermediates)
            mrich.var("products", products)
            mrich.var("reactions", reactions)

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
            self._product_poses._name = f"product poses of {self}"
        return self._product_poses

    @property
    def product_compounds(self) -> "CompoundSet":
        if self._product_compounds is None:
            self._product_compounds = self.products.compounds
            self._product_compounds._name = f"products of {self}"
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
                    price=str(ingredient.price),
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
                price=str(ingredient.price),
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
                mrich.error(f"problem w/ node {key=}")
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
                mrich.error(f"problem w/ edge {s=} {t=}")
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
            try:
                title = f"Recipe<br><sup>price={self.price}</sup>"
            except AssertionError:
                title = f"Recipe"

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

        mrich.header("Recipe")

        if price:
            price = self.price
            if price:
                print("\nprice", price.amount, dict(unit=price.currency))
                # mrich.var('lead-time', self.lead_time, dict(unit='working days'))

        mrich.var("\n#products", len(self.products))

        if len(self.products) < 100:
            for product in self.products:
                print(str(product.compound), f"{product.amount:.2f}", dict(unit="mg"))

        mrich.var("\n#intermediates", len(self.intermediates))
        if len(self.intermediates) < 100:
            for intermediate in self.intermediates:
                print(
                    str(intermediate.compound),
                    f"{intermediate.amount:.2f}",
                    dict(unit="mg"),
                )

        mrich.var("\n#reactants", len(self.reactants))
        if len(self.reactants) < 100:
            for reactant in self.reactants:
                print(str(reactant.compound), f"{reactant.amount:.2f}", dict(unit="mg"))

        mrich.var("\n#reactions", len(self.reactions))
        if len(self.reactions) < 100:
            for reaction in self.reactions:
                print(str(reaction), reaction.reaction_str, dict(unit=reaction.type))

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

    def write_json(self, file, *, extra: dict | None = None, indent="\t", **kwargs):
        """Serialise this recipe object and write it to disk

        :param file:
        :param extra:
        :param indent:  (Default value = '\t')

        """
        import json

        data = self.get_dict(serialise_price=True, **kwargs)

        if extra:
            data.update(extra)

        mrich.writing(file)
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
        try:
            if price and serialise_price:
                data["price"] = self.price.get_dict()
            elif price:
                data["price"] = self.price
        except AssertionError as e:
            mrich.warning(f"Could not get price: {e}")
            data["price"] = None

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

    def get_routes(self):
        return self.products.get_routes(permitted_reactions=self.reactions)

    def write_CAR_csv(
        self, file: "str | Path", return_df: bool = False
    ) -> "DataFrame | None":
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
        from pathlib import Path

        # solve each product's reaction

        file = str(Path(file).resolve())

        rows = []

        routes = self.get_routes()

        for sub_recipe in routes:

            product = sub_recipe.product

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
            mrich.writing(this_file)
            subset.to_csv(this_file, index=False)

        mrich.writing(file)
        df.to_csv(file, index=False)

        return df

    def write_reactant_csv(
        self, file: "str | Path", return_df: bool = False
    ) -> "DataFrame | None":
        """Detailed CSV output including reactant information for purchasing and information on the downstream synthetic use

        Reactant
        ========

        - ID
        - SMILES
        - Inchikey

        Quote
        =====

        - Supplier
        - Catalogue
        - Entry
        - Lead-time
        - Quoted amount
        - Quote currency
        - Quote price
        - Quote purity

        Downstream
        ==========

        - num_reaction_dependencies
        - num_product_dependencies
        - reaction_dependencies
        - product_dependencies

        """
        # - remove_with

        from pandas import DataFrame

        # from rich import print
        from .cset import CompoundSet
        from .rset import ReactionSet

        data = []

        routes = self.get_routes()

        for reactant in mrich.track(
            self.reactants, prefix="Constructing reactant DataFrame"
        ):
            quote = reactant.quote

            d = dict(
                hippo_id=reactant.compound_id,
                smiles=reactant.smiles,
                inchikey=reactant.inchikey,
                required_amount_mg=reactant.amount,
                quoted_amount=quote.amount,
                quote_currency=quote.currency,
                quote_price=quote.price.amount,
                quote_lead_time_days=quote.lead_time,
                quote_supplier=quote.supplier,
                quote_catalogue=quote.catalogue,
                quote_entry=quote.entry,
                quoted_smiles=quote.smiles,
                quoted_purity=quote.purity,
            )

            downstream_routes = []
            downstream_reactions = []

            for route in routes:
                if reactant in route.reactants:
                    downstream_routes.append(route)
                for reaction in route.reactions:
                    if reactant in reaction.reactants:
                        downstream_reactions.append(reaction)

            downstream_products = CompoundSet(
                self.db, set(route.product.id for route in downstream_routes)
            )
            downstream_reactions = ReactionSet(
                self.db, set(reaction.id for reaction in downstream_reactions)
            )

            def get_scaffold_series():

                bases = downstream_products.bases

                if not bases:
                    bases = downstream_products[0:]

                return bases.ids

            d["num_reaction_dependencies"] = len(downstream_reactions)
            d["num_product_dependencies"] = len(downstream_products)
            d["reaction_dependencies"] = downstream_reactions.ids
            d["product_dependencies"] = downstream_products.ids
            d["chemistry_types"] = ", ".join(set(downstream_reactions.types))
            d["scaffold_series"] = get_scaffold_series()

            data.append(d)

        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df

        return None

    def write_product_csv(self, file: "str | Path", return_df: bool = False):
        """Detailed CSV output including product information for selection and synthesis"""

        from pandas import DataFrame

        # from rich import print
        from .cset import CompoundSet
        from .rset import ReactionSet

        data = []

        routes = self.get_routes()

        pose_map = self.db.get_compound_id_pose_ids_dict(self.products.compounds)

        for product in mrich.track(
            self.products, prefix="Constructing product DataFrame"
        ):

            d = dict(
                hippo_id=product.compound_id,
                smiles=product.smiles,
                inchikey=product.inchikey,
                required_amount_mg=product.amount,
            )

            upstream_routes = []
            upstream_reactions = []

            for route in routes:
                if product in route.products:
                    upstream_routes.append(route)

                    for reaction in route.reactions:
                        upstream_reactions.append(reaction)

            upstream_reactions = ReactionSet(
                self.db, set(reaction.id for reaction in upstream_reactions)
            )

            def get_scaffold_series():

                if bases := product.bases:
                    return bases.ids, False

                else:
                    return [product.id], True

            poses = pose_map.get(product.id, set())

            d["num_poses"] = len(poses)
            d["poses"] = poses
            d["tags"] = product.tags
            d["num_routes"] = len(upstream_routes)
            d["num_reaction_steps"] = set(
                len(route.reactions) for route in upstream_routes
            )
            d["reaction_dependencies"] = upstream_reactions.ids
            d["reactant_dependencies"] = set(
                sum([route.reactants.ids for route in upstream_routes], [])
            )
            d["route_ids"] = [route.id for route in upstream_routes]
            d["chemistry_types"] = ", ".join(set(upstream_reactions.types))
            series, is_base = get_scaffold_series()
            d["is_scaffold"] = is_base
            d["scaffold_series"] = series

            data.append(d)

        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df

        return None

    def write_chemistry_csv(self, file: "str | Path", return_df: bool = True):
        """Detailed CSV output synthetis information for chemistry types in this set"""

        from pandas import DataFrame

        from rich import print
        from .cset import CompoundSet
        from .rset import ReactionSet

        data = []

        # get compounds

        scaffolds = CompoundSet(self.db)

        for product in self.products:

            if bases := product.bases:
                scaffolds += bases
            else:
                scaffolds.add(product.compound)

        routes = self.get_routes()

        route_types = {}

        for compound in scaffolds:

            elabs = (
                self.products.compounds.get_by_base(base=compound, none="quiet") or []
            )

            d = dict(
                scaffold_id=compound.id,
                product_id=compound.id,
                smiles=compound.smiles,
                inchikey=compound.inchikey,
                num_elaborations=len(elabs),
                is_scaffold=True,
            )

            upstream_routes = []
            for route in routes:
                if compound in route.products:
                    upstream_routes.append(route)

            if not upstream_routes:
                mrich.warning(f"No routes to scaffold={compound}")
                continue

            d["num_routes"] = len(upstream_routes)

            for j, route in enumerate(upstream_routes):
                d[f"route_{j+1}_num_steps"] = len(route.reactions)

                group = route_types.setdefault(compound.id, set())
                group.add(tuple([r.type for r in route.reactions]))

                for k, reaction in enumerate(route.reactions):
                    key = f"route_{j+1}_reaction_{k+1}"

                    product = reaction.product

                    d[f"{key}_type"] = reaction.type
                    d[f"{key}_product_smiles"] = product.smiles
                    d[f"{key}_product_id"] = product.id
                    d[f"{key}_product_yield"] = reaction.product_yield

                    for i, reactant in enumerate(reaction.reactants):
                        d[f"{key}_reactant_{i+1}_smiles"] = reactant.smiles
                        d[f"{key}_reactant_{i+1}_id"] = reactant.id

            data.append(d)

        missing_bases = {}

        for compound in self.products.compounds:

            if compound in scaffolds:
                continue

            upstream_routes = []
            for route in routes:
                if compound in route.products:
                    upstream_routes.append(route)

            bases = compound.bases

            for base in bases:

                if base.id not in route_types:
                    group = missing_bases.setdefault(base.id, [])
                    group.append(compound.id)
                    continue

                else:
                    for route in upstream_routes:
                        chem_types = tuple([r.type for r in route.reactions])

                        if chem_types not in route_types[base.id]:
                            mrich.success(base)
                            mrich.success(chem_types)
                            raise ValueError(
                                "Scaffold has route not present in dataframe"
                            )

        for base_id, elab_ids in missing_bases.items():

            compound = self.db.get_compound(id=sorted(elab_ids)[0])

            d = dict(
                scaffold_id=base_id,
                product_id=compound.id,
                smiles=compound.smiles,
                inchikey=compound.inchikey,
                num_elaborations=len(elab_ids),
                is_scaffold=False,
            )

            upstream_routes = []
            for route in routes:
                if compound in route.products:
                    upstream_routes.append(route)

            if not upstream_routes:
                mrich.error(f"No routes to elab {compound}")
                raise ValueError(f"No routes to elab {compound}")

            d["num_routes"] = len(upstream_routes)

            for j, route in enumerate(upstream_routes):
                d[f"route_{j+1}_num_steps"] = len(route.reactions)

                group = route_types.setdefault(compound.id, set())
                group.add(tuple([r.type for r in route.reactions]))

                for k, reaction in enumerate(route.reactions):
                    key = f"route_{j+1}_reaction_{k+1}"

                    product = reaction.product

                    d[f"{key}_type"] = reaction.type
                    d[f"{key}_product_smiles"] = product.smiles
                    d[f"{key}_product_id"] = product.id
                    d[f"{key}_product_yield"] = reaction.product_yield

                    for i, reactant in enumerate(reaction.reactants):
                        d[f"{key}_reactant_{i+1}_smiles"] = reactant.smiles
                        d[f"{key}_reactant_{i+1}_id"] = reactant.id

            data.append(d)

        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df

        return None

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
            mrich.warning(f"{null_count} poses have no fingerprint")

        return poses.fingerprint

    def __flag_modification(self):
        self._product_interactions = None
        self._score = None
        self._product_compounds = None
        self._product_poses = None

    def check_integrity(self, debug: bool = False) -> bool:

        reaction_products = self.reactions.products
        reaction_reactants = self.reactions.reactants

        # all products should have a reaction
        for product in self.products:
            if product not in reaction_products:
                mrich.error(f"Product: {product} does not have associated reaction")
                return False

        # all intermediates should be the product of a reaction and the reactant of a reaction
        for intermediate in self.intermediates:
            if intermediate not in reaction_products:
                mrich.error(
                    f"Intermediate: {intermediate} is not the product of any included reactions"
                )
                return False
            if intermediate not in reaction_reactants:
                mrich.error(
                    f"Intermediate: {intermediate} is not a reactant of any included reactions"
                )
                return False

        # all reactions should have enough reactant

        for reaction in self.reactions:

            product_ingredient = self.products(compound_id=reaction.product_id)

            if product_ingredient is None:
                product_ingredient = self.intermediates(compound_id=reaction.product_id)

            if debug and reaction.product_yield < 1.0:
                mrich.debug(f"{reaction}.product_yield={reaction.product_yield}")

            for reactant in reaction.reactants:

                reactant_ingredient = self.intermediates(compound_id=reactant.id)

                if reactant_ingredient is None:
                    reactant_ingredient = self.reactants(compound_id=reactant.id)

                required_amount = product_ingredient.amount / reaction.product_yield

                if reactant_ingredient.amount < required_amount:
                    mrich.error(
                        f"Not enough of {reactant_ingredient.compound}: {reactant_ingredient.amount} < {required_amount}"
                    )
                    return False

        return True

    ### DUNDERS

    def __str__(self):

        if self.score:
            s = f"score={self.score:.3f}"
        else:
            s = ""

        if self.hash:
            return f"Recipe_{self.hash}({s})"

        return f"Recipe({s})"

    def __longstr(self) -> str:
        """Unformatted string representation"""
        if self.intermediates:
            s = f"{self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions}"
        else:
            s = f"{self.reactants} --> {self.products} via {self.reactions}"

        if self.score:
            s += f", score={self.score:.3f}"

        if self.hash:
            return f"Recipe_{self.hash}({s})"

        return f"Recipe({s})"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.__longstr()}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self.__longstr()}"

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

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"Route #{self.id}: {self.product_compound}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"


class RouteSet:
    """A set of Route objects"""

    def __init__(self, db, routes):

        data = {}
        for route in routes:
            # assert isinstance(route, Route)
            data[route.id] = route

        self._data = data
        self._db = db
        self._cluster_map = None
        self._permitted_clusters = None
        self._current_cluster = None

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
        for d in mrich.track(data["routes"].values(), prefix="Loading Routes..."):
            route_id = d["id"]
            new_data[route_id] = Route.from_json(db=db, path=None, data=d)

        self._data = new_data
        self._db = db
        self._cluster_map = None
        self._permitted_clusters = None
        self._current_cluster = None

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
    def products(self) -> "CompoundSet":
        """Return a :class:`.CompoundSet` of all the route products"""
        from .cset import CompoundSet

        return CompoundSet(self.db, self.product_ids)

    @property
    def str_ids(self):
        """Return an SQL formatted tuple string of the :class:`.Route` ID's"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def ids(self):
        """Return the :class:`.Route` IDs"""
        return self.data.keys()

    @property
    def cluster_map(self) -> dict[tuple, set]:
        """Create a dictionary grouping routes by their scaffold/base cluster.

        :returns: A dictionary mapping a tuple of scaffold :class:`.Compound` IDs to a set of :class:`.Route` ID's to their superstructures.
        """

        if self._cluster_map is None:

            # get route mapping
            pairs = self.db.select_where(
                query="route_product, route_id",
                key=f"route_id IN {self.str_ids}",
                table="route",
                multiple=True,
            )

            route_map = {route_product: route_id for route_product, route_id in pairs}

            # group compounds by cluster
            compound_clusters = self.db.get_compound_cluster_dict(cset=self.products)

            # create the map
            self._cluster_map = {}
            for cluster, compounds in compound_clusters.items():
                self._cluster_map[cluster] = []
                for compound in compounds:
                    route_id = route_map.get(compound, None)
                    if not route_id:
                        continue
                    self._cluster_map[cluster].append(route_id)

                if not self._cluster_map[cluster]:
                    del self._cluster_map[cluster]

        return self._cluster_map

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

    def balanced_pop(
        self, permitted_clusters: set[tuple] | None = None, debug: bool = False
    ) -> "Route":

        if not self._data:
            mrich.print("RouteSet depleted")
            return None

        if not self.cluster_map:
            # mrich.warning("RouteSet.cluster_map depleted but _data isn't...")
            return self.pop()

        # store the permitted clusters (or all clusters) list as property

        if self._permitted_clusters is None:
            if permitted_clusters:
                permitted_clusters = set(
                    (cluster,) if isinstance(cluster, int) else cluster
                    for cluster in permitted_clusters
                )

                self._permitted_clusters = []
                for cluster in permitted_clusters:
                    if cluster not in self.cluster_map:
                        mrich.warning(
                            cluster, "in permitted_clusters but not cluster_map"
                        )
                    else:
                        self._permitted_clusters.append(cluster)

            else:
                self._permitted_clusters = list(self.cluster_map.keys())

        if self._current_cluster is None:
            self._current_cluster = self._permitted_clusters[0]

        ### pop a Route

        if debug:
            mrich.debug(f"Would pop Route from {self._current_cluster=}")

        cluster = self._current_cluster

        # pop the last route id from the given cluster

        try:
            route_id = self.cluster_map[cluster].pop()
        except IndexError:
            mrich.print(self._permitted_clusters)
            mrich.print(self.cluster_map)
            raise
        except AttributeError:
            mrich.print(cluster)
            mrich.print(self.cluster_map)
            raise
        except KeyError:
            mrich.print("cluster", cluster)
            mrich.print("self._permitted_clusters", self._permitted_clusters)
            mrich.print("self.cluster_map.keys()", self.cluster_map.keys())
            raise

        # clean up empty clusters

        if debug:
            mrich.debug("Popped route", route_id)

        # get the Route object

        if route_id in self._data:
            route = self._data[route_id]
            del self._data[route_id]
        else:
            # if debug:
            mrich.debug("Route not present")
            return self.balanced_pop()

        ### increment cluster

        # def increment_cluster(cluster):
        n = len(self._permitted_clusters)
        if n > 1:
            for i, cluster in enumerate(self._permitted_clusters):
                if cluster == self._current_cluster:
                    if i == n - 1:
                        self._current_cluster = self._permitted_clusters[0]
                    else:
                        self._current_cluster = self._permitted_clusters[i + 1]
                    break
            else:
                raise IndexError("This should never be reached...")

        # increment_cluster()

        if not self.cluster_map[cluster]:
            del self.cluster_map[cluster]
            if not self.cluster_map:
                mrich.debug("RouteSet.cluster_map depleted")
            self._permitted_clusters = [
                c for c in self._permitted_clusters if c != cluster
            ]
            # if debug:
            mrich.debug("Depleted cluster", cluster)

            if not self._permitted_clusters:
                mrich.debug("Depleted all permitted clusters", cluster)
                mrich.debug("Removing cluster restriction", cluster)
                self._permitted_clusters = list(self.cluster_map.keys())
                self._current_cluster = None

        if debug:
            mrich.debug("#Routes in set", len(self._data))

        return route

    def shuffle(self):
        """Randomly shuffle the routes in this set"""
        import random

        items = list(self.data.items())
        random.shuffle(items)
        self._data = dict(items)

        ### shuffle the cluster map as well

        for cluster, routes in self.cluster_map.items():
            random.shuffle(routes)
            self.cluster_map[cluster] = routes

    ### DUNDERS

    def __len__(self):
        return len(self.data)

    def __str__(self) -> str:
        """Unformatted string representation"""
        return "{" f"Route × {len(self)}" "}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"

    def __iter__(self):
        return iter(self.data.values())


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
        from json import JSONDecodeError

        self._db = db
        self._json_directory = Path(directory)
        self._json_pattern = pattern

        self._json_paths = {}
        for path in self._json_directory.glob(self._json_pattern):
            self._json_paths[
                path.name.removeprefix("Recipe_").removesuffix(".json")
            ] = path.resolve()

        mrich.reading(f"{directory}/{pattern}")

        self._recipes = {}
        for key, path in mrich.track(
            self._json_paths.items(), prefix="Loading recipes"
        ):
            try:
                recipe = Recipe.from_json(
                    db=self.db,
                    path=path,
                    allow_db_mismatch=True,
                    debug=False,
                    db_mismatch_warning=False,
                )
            except JSONDecodeError:
                mrich.error(f"Bad JSON in {path}")
                continue
            recipe._hash = key
            self._recipes[key] = recipe

        mrich.success("Loaded", len(self), "Recipes")

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
            recipes = mrich.track(recipes, prefix=f"Calculating {self} values...")

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

    def keys(self):
        return self._recipes.keys()

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
                mrich.error(
                    f"Unsupported type for RecipeSet.__getitem__(): {key=} {type(key)}"
                )

        return None

    def __iter__(self):
        return iter(self._recipes.values())

    def __contains__(self, key):
        assert isinstance(key, str)
        return key in self._recipes

    def __str__(self) -> str:
        """Unformatted string representation"""
        return "{" f"Recipe × {len(self)}" "}"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"
