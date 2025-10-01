from dataclasses import dataclass, field

from .compound import Ingredient

import mcol

import mrich


class Recipe:
    """A Recipe stores data corresponding to a specific synthetic recipe involving several products, reactants, intermediates, and reactions."""

    _db = None

    def __init__(
        self,
        db: "Database",
        *,
        products: "IngredientSet | None" = None,
        reactants: "IngredientSet | None" = None,
        intermediates: "IngredientSet | None" = None,
        reactions: "ReactionSet | None" = None,
        compounds: "IngredientSet | None" = None,
    ):

        from .cset import IngredientSet
        from .rset import ReactionSet

        if products is None:
            products = IngredientSet(db)

        if reactants is None:
            reactants = IngredientSet(db)

        if intermediates is None:
            intermediates = IngredientSet(db)

        if compounds is None:
            compounds = IngredientSet(db)

        if reactions is None:
            reactions = ReactionSet(db)

        # check typing
        assert isinstance(products, IngredientSet)
        assert isinstance(reactants, IngredientSet)
        assert isinstance(intermediates, IngredientSet)
        assert isinstance(compounds, IngredientSet)
        assert isinstance(reactions, ReactionSet)

        self._products = products
        self._reactants = reactants
        self._intermediates = intermediates
        self._reactions = reactions
        self._compounds = compounds
        self._db = db
        self._hash = None

        self._score = None

        # caches
        self._product_compounds = None
        self._poses = None
        self._interactions = None
        self._combined_compounds = None

    ### FACTORIES

    @classmethod
    def from_reaction(
        cls,
        reaction,
        amount=1,
        *,
        debug: bool = False,
        pick_cheapest: bool = True,
        permitted_reactions: "ReactionSet | None" = None,
        quoted_only: bool = False,
        supplier: None | str = None,
        unavailable_reaction: str = "error",
        reaction_checking_cache: dict[int, bool] = None,
        reaction_reactant_cache: dict[int, bool] = None,
        inner: bool = False,
        get_ingredient_quotes: bool = True,
    ) -> "Recipe | list[Recipe]":
        """Create a :class:`.Recipe` from a :class:`.Reaction` and its upstream dependencies

        :param reaction: reaction to create recipe from
        :param amount: amount in ``mg`` (Default value = 1)
        :param debug: bool: increase verbosity for debugging (Default value = False)
        :param pick_cheapest: bool: choose the cheapest solution (Default value = True)
        :param permitted_reactions: once consider reactions in this set (Default value = None)
        :param quoted_only: bool: only allow reactants with quotes (Default value = False)
        :param supplier: None | str: optionally restrict quotes to only this supplier (Default value = None)
        :param unavailable_reaction: define the behaviour for when a reaction has unavailable reactants (Default value = 'error')
        :param inner: used to indicate that this is a recursive call (Default value = False)
        :param get_ingredient_quotes: get quotes for ingredients in this recipe

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
        reactions: "ReactionSet",
        amount: float = 1,
        pick_cheapest: bool = True,
        permitted_reactions: "ReactionSet | None" = None,
        final_products_only: bool = True,
        return_products: bool = False,
        supplier: str | None = None,
        use_routes: bool = False,
        debug: bool = False,
        **kwargs,
    ) -> "Recipe | list[Recipe] | CompoundSet":
        """Create a :class:`.Recipe` from a :class:`.ReactionSet` and its upstream dependencies

        :param reactions: reactions to create recipe from
        :param amount: amount in ``mg`` (Default value = 1)
        :param debug: bool: increase verbosity for debugging (Default value = False)
        :param pick_cheapest: bool: choose the cheapest solution (Default value = True)
        :param permitted_reactions: once consider reactions in this set (Default value = None)
        :param final_products_only: don't get routes to intermediates (Default value = True)
        :param return_products: return the :class:`.CompoundSet` of products instead (Default value = False)

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
            supplier=supplier,
            use_routes=use_routes,
            **kwargs,
        )

        return recipe

    @classmethod
    def from_compounds(
        cls,
        compounds: "CompoundSet",
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
        unavailable_reaction: str = "error",
        reaction_checking_cache: dict[int, bool] | None = None,
        reaction_reactant_cache: dict[int, bool] | None = None,
        use_routes: bool = False,
        **kwargs,
    ):
        """Create recipe(s) to synthesis products in the :class:`.CompoundSet`

        :param compounds: set of compounds to find routes for
        :param solve_combinations: bool: combinatorially combine all individual routes (Default value = True)
        :param pick_first: return the first solution without comparison (Default value = False)
        :param warn_multiple_solutions: warn if a compound has multiple routes (Default value = True)
        :param pick_cheapest_inner_routes: for each compound choose the cheapest route (Default value = False)
        :param reaction: reaction to create recipe from
        :param amount: amount in ``mg`` (Default value = 1)
        :param debug: bool: increase verbosity for debugging (Default value = False)
        :param pick_cheapest: bool: choose the cheapest solution (Default value = True)
        :param permitted_reactions: once consider reactions in this set (Default value = None)
        :param quoted_only: bool: only allow reactants with quotes (Default value = False)
        :param supplier: None | str: optionally restrict quotes to only this supplier (Default value = None)
        :param unavailable_reaction: define the behaviour for when a reaction has unavailable reactants (Default value = 'error')

        """

        from .cset import CompoundSet

        assert isinstance(compounds, CompoundSet)

        # if permitted_reactions:
        #   raise NotImplementedError

        db = compounds.db

        n_comps = len(compounds)

        assert n_comps

        if not hasattr(amount, "__iter__"):
            amount = [amount] * n_comps

        if use_routes:
            route_lookup = db.get_product_id_routes_dict()

            if supplier:
                raise NotImplementedError
                # supplier_lookup = db.get_compound_id_suppliers_dict()

        options = []

        ok = 0
        mrich.var("#compounds", n_comps)

        for comp, a in mrich.track(
            zip(compounds, amount),
            prefix="Solving individual compound recipes...",
            total=n_comps,
        ):
            comp_options = []

            if use_routes:

                if comp.id not in route_lookup:
                    mrich.error("No routes to", comp)
                    continue

                comp_options = []
                for route_id in route_lookup[comp.id]:
                    route = db.get_route(id=route_id)
                    comp_options.append(route)

            else:

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
                    mrich.error(
                        f"No solutions for compound={comp} ({comp.reactions.ids=})"
                    )
                    continue

            if pick_cheapest and len(comp_options) > 1:
                if warn_multiple_solutions:
                    mrich.warning(
                        f"Multiple solutions for", comp, "(", len(comp_options), ")"
                    )
                if debug:
                    mrich.debug("Picking cheapest...")
                priced = [r for r in comp_options if r.price]
                comp_options = sorted(priced, key=lambda r: r.price)[:1]

            if warn_multiple_solutions and len(comp_options) > 1:
                mrich.warning(f"Multiple solutions for compound={comp}")
                if debug:
                    mrich.debug(f"{comp_options=}")
            else:
                if n_comps <= 200:
                    mrich.success(f"Found solution for compound={comp}")
                ok += 1
                mrich.set_progress_field("ok", ok)
                mrich.set_progress_field("n", n_comps)

            options.append(comp_options)

        assert all(options)

        from itertools import product

        mrich.print("Solving recipe combinations...")
        combinations = list(product(*options))

        if not solve_combinations:
            return combinations

        # if pick_first:
        #     combinations = [combinations[0]]

        solutions = []

        if n_comps > 1:
            generator = mrich.track(
                combinations, prefix="Combining recipes...", total=len(combinations)
            )
        else:
            generator = combinations

        ok = 0
        for combo in generator:

            if debug:
                mrich.debug(f"Combination of {len(combo)} recipes")

            if not combo:
                continue

            solution = combo[0]

            for i, recipe in enumerate(combo[1:]):
                if debug:
                    mrich.debug(i + 1)
                solution += recipe

            solutions.append(solution)
            ok += 1
            mrich.set_progress_field("ok", ok)
            mrich.set_progress_field("n", len(combinations))

        if not solutions:
            mrich.error("No solutions")
            return None

        if pick_first:
            return solutions[0]

        if pick_cheapest:
            mrich.debug("Calculating prices...")
            priced = [r for r in solutions if r.price]
            mrich.print("Picking cheapest from", len(priced), "options")
            if not priced:
                mrich.error("0 recipes with prices, can't choose cheapest")
                return solutions
            return sorted(priced, key=lambda r: r.price)[0]

        return solutions

    @classmethod
    def from_reactants(
        cls,
        reactants: "CompoundSet | IngredientSet",
        amount: float = 1,
        debug: bool = False,
        return_products: bool = False,
        supplier: str | None = None,
        pick_cheapest: bool = False,
        use_routes: bool = False,
        **kwargs,
    ) -> "list[Recipe] | Recipe | CompoundSet":
        """Find the maximal recipe from a given set of reactants

        :param reactants: :class:`.CompoundSet` or :class:`.IngredientSet` for the reactants. Ingredient amounts are ignored
        :param amount: amount of each product needed (Default value = 1)
        :param debug: increase verbosity (Default value = False)
        :param return_products: return products instead of recipe (Default value = False)
        :param kwargs: passed to :meth:`.Recipe.from_reactions`

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
            debug=debug,
            return_products=return_products,
            supplier=supplier,
            use_routes=use_routes,
            **kwargs,
        )

        return recipe

    @classmethod
    def from_json(
        cls,
        db: "Database",
        path: "str | Path",
        debug: bool = True,
        allow_db_mismatch: bool = False,
        clear_quotes: bool = False,
        data: dict = None,
        db_mismatch_warning: bool = True,
    ):
        """Load a serialised recipe from a JSON file

        :param db: database to link
        :param path: path to JSON
        :param debug: increase verbosity (Default value = True)
        :param allow_db_mismatch: allow a database mismatch (Default value = False)
        :param clear_quotes: ignore reactant quotes (Default value = False)
        :param data: serialised data (Default value = None)

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

        if "compounds" in data:
            compounds = IngredientSet.from_ingredient_dicts(
                db, data["compounds"], supplier=data["compound_supplier"]
            )
        else:
            compounds = IngredientSet(db)

        if clear_quotes:
            reactants.df["quote_id"] = None
            reactants.df["quoted_amount"] = None
            compounds.df["quote_id"] = None
            compounds.df["quoted_amount"] = None

        # ReactionSet
        reactions = ReactionSet(db, data["reaction_ids"], sort=False)

        if debug:
            mrich.var("reactants", reactants)
            mrich.var("intermediates", intermediates)
            mrich.var("products", products)
            mrich.var("reactions", reactions)
            mrich.var("compounds", compounds)

        # Create the object
        self = cls.__new__(cls)
        self.__init__(
            db,
            products=products,
            reactants=reactants,
            intermediates=intermediates,
            reactions=reactions,
            compounds=compounds,
        )

        return self

    ### PROPERTIES

    @property
    def db(self) -> "Database":
        """Associated :class:`.Database:"""
        return self._db

    @property
    def products(self) -> "IngredientSet":
        """Product :class:`.IngredientSet`"""
        return self._products

    @property
    def compounds(self) -> "IngredientSet":
        """Product :class:`.IngredientSet`"""
        return self._compounds

    @compounds.setter
    def compounds(self, a: "IngredientSet"):
        """Set the compounds"""
        self._compounds = a
        self.__flag_modification()

    @property
    def poses(self) -> "PoseSet":
        """Product poses"""
        if self._poses is None:
            self._poses = self.combined_compounds.poses
            self._poses._name = f"poses of {self}"
        return self._poses

    @property
    def product_compounds(self) -> "CompoundSet":
        """Product compounds"""
        if self._product_compounds is None:
            self._product_compounds = self.products.compounds
            self._product_compounds._name = f"products of {self}"
        return self._product_compounds

    @property
    def combined_compound_ids(self) -> set[int]:
        return set(self.product_compounds.ids) | set(self.compounds.ids)

    @property
    def combined_compounds(self) -> "CompoundSet":
        """Combined product and no-chem compounds"""
        if self._combined_compounds is None:
            from .cset import CompoundSet

            self._combined_compounds = CompoundSet(self.db, self.combined_compound_ids)
            self._combined_compounds._name = f"combined compounds of {self}"
        return self._combined_compounds

    @property
    def interactions(self) -> "InteractionSet":
        """Product pose interactions"""
        if self._interactions is None:
            self._interactions = self.poses.interactions
        return self._interactions

    @property
    def product(self) -> "Ingredient":
        """Return single product (if there's only one)"""
        assert len(self.products) == 1
        return self.products[0]

    @products.setter
    def products(self, a: "IngredientSet"):
        """Set the products"""
        self._products = a
        self.__flag_modification()

    @property
    def reactants(self):
        """Reactant :class:`.IngredientSet`"""
        return self._reactants

    @reactants.setter
    def reactants(self, a: "IngredientSet"):
        """Set the reactants"""
        self._reactants = a
        self.__flag_modification()

    @property
    def intermediates(self) -> "IngredientSet":
        """Intermediates :class:`.IngredientSet`"""
        return self._intermediates

    @intermediates.setter
    def intermediates(self, a: "IngredientSet"):
        """Set the intermediates"""
        self._intermediates = a
        self.__flag_modification()

    @property
    def reactions(self) -> "ReactionSet":
        """Intermediates :class:`.IngredientSet`"""
        return self._reactions

    @reactions.setter
    def reactions(self, a: "ReactionSet"):
        """Set the reactions"""
        self._reactions = a
        self.__flag_modification()

    @property
    def price(self) -> "Price":
        """Get the price of the reactants"""
        return self.reactants.get_price() + self.compounds.get_price()

    @property
    def num_products(self) -> int:
        """Return the number of products"""
        return len(self.products)

    @property
    def num_compounds(self) -> int:
        """Return the number of compounds"""
        return len(self.combined_compound_ids)

    @property
    def num_reactions(self):
        """Return the number of reactions"""
        return len(self.reactions)

    @property
    def num_reactants(self):
        """Return the number of reactants"""
        return len(self.reactants)

    @property
    def num_intermediates(self):
        """Return the number of intermediates"""
        return len(self.intermediates)

    @property
    def hash(self) -> str:
        """Return the unique hash string"""
        return self._hash

    @property
    def score(self):
        """Return the Recipe score"""
        return self._score

    @property
    def type(self) -> str:

        if self.empty:
            return "EMPTY"

        chem = bool(self.reactions)
        nochem = bool(self.compounds)

        if chem and nochem:
            return "MIXED"

        if chem and not nochem:
            return "CHEM"

        if nochem and not chem:
            return "NOCHEM"

    @property
    def empty(self) -> bool:
        """Is this Recipe empty?"""

        if self.reactants:
            return False

        if self.products:
            return False

        if self.intermediates:
            return False

        if self.reactions:
            return False

        if self.compounds:
            return False

        return True

    ### METHODS

    def get_price(self, supplier: str | None = None) -> "Price":
        """get the reactants price. See :meth:`.IngredientSet.get_price`

        :param supplier: restrict quotes to this supplier

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

    def sankey(self, title: str | None = None) -> "graph_objects.Figure":
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

    def summary(self, price: bool = True) -> None:
        """Print a summary of this recipe

        :param price: print the price (Default value = True)

        """

        import mcol

        mrich.h1(str(self))

        if price:
            price = self.price
            if price:
                mrich.var("\nprice", price.amount, price.currency)
                # mrich.var('lead-time', self.lead_time, 'working days))

        if self.products:
            mrich.h3(f"{len(self.products)} products")

            if len(self.products) < 100:
                for product in self.products:
                    mrich.var(str(product.compound), f"{product.amount:.2f}", "mg")

        if self.intermediates:
            mrich.h3(f"{len(self.intermediates)} intermediates")

            if len(self.intermediates) < 100:
                for intermediate in self.intermediates:
                    mrich.var(
                        str(intermediate.compound),
                        f"{intermediate.amount:.2f}",
                        "mg",
                    )

        if self.reactants:
            mrich.h3(f"{len(self.reactants)} reactants")

            if len(self.reactants) < 100:
                for reactant in self.reactants:
                    mrich.var(str(reactant.compound), f"{reactant.amount:.2f}", "mg")

        if self.reactions:
            mrich.h3(f"{len(self.reactions)} reactions")

            if len(self.reactions) < 100:
                for reaction in self.reactions:
                    mrich.var(str(reaction), reaction.reaction_str, reaction.type)

        if self.compounds:

            mrich.h3(f"{len(self.compounds)} compounds")

            if len(self.compounds) < 100:
                for compound in self.compounds:
                    mrich.var(str(compound.compound), f"{compound.amount:.2f}", "mg")

    def get_ingredient(self, id) -> "Ingredient":
        """Get an ingredient by its compound ID

        :param id: compound ID

        """
        matches = [r for r in self.reactants if r.id == id]
        if not matches:
            matches = [r for r in self.intermediates if r.id == id]
        if not matches:
            matches = [r for r in self.products if r.id == id]

        assert len(matches) == 1
        return matches[0]

    def add_to_all_reactants(self, amount: float = 20) -> None:
        """Increment all reactants by this amount

        :param amount: amount in ``mg`` (Default value = 20)

        """
        self.reactants.df["amount"] += amount

    def write_json(
        self,
        file: "str | Path",
        *,
        extra: dict | None = None,
        indent: str = "\t",
        **kwargs,
    ) -> None:
        """Serialise this recipe object and write it to disk

        :param file: write to this path
        :param extra: extra data to serialise
        :param indent: indentation whitespace (Default value = '\t')

        """
        import json
        from pathlib import Path

        file = Path(file).resolve()

        assert file.parent.exists()

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
        compound_supplier: bool = True,
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

        :param price: include the price (Default value = True)
        :param reactant_supplier: include the supplier (Default value = True)
        :param database: include the database (Default value = True)
        :param timestamp: add a timestamp (Default value = True)
        :param compound_ids_only: ID's only (instead of full :attr:`.IngredientSet.df`) (Default value = False)
        :param products: include products (Default value = True)
        :param serialise_price: serialise :class:`.Price` object (Default value = False)

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

        if compound_supplier:
            data["compound_supplier"] = self.compounds.supplier

        # IngredientSets
        if compound_ids_only:
            data["reactant_ids"] = self.reactants.compound_ids
            data["intermediate_ids"] = self.intermediates.compound_ids
            if products:
                data["products_ids"] = self.products.compound_ids
            data["compound_ids"] = self.compounds.compound_ids

        else:
            data["reactants"] = self.reactants.df.to_dict(orient="list")
            data["intermediates"] = self.intermediates.df.to_dict(orient="list")
            if products:
                data["products"] = self.products.df.to_dict(orient="list")
            data["compounds"] = self.compounds.df.to_dict(orient="list")

        # ReactionSet
        data["reaction_ids"] = self.reactions.ids

        return data

    def get_routes(self, return_ids: bool = False) -> "RouteSet":
        """Get routes"""
        return self.products.get_routes(
            permitted_reactions=self.reactions, return_ids=return_ids
        )

    def write_CAR_csv(
        self, file: "str | Path", return_df: bool = False
    ) -> "DataFrame | None":
        """Prepares CSVs for use with CAR.

        .. attention::

            This method requires a populated `route` table. For a workaround use :meth:`.CompoundSet.write_CAR_csv` instead

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

        :param file: file to write to
        :param return_df: return the dataframe (Default value = False)

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
                        # mrich.warning(f"More than two reactants for {reaction=}")
                        for j, r in enumerate(reaction.reactants):
                            row[f"reactant-{j+1}-{i}"] = reaction.reactants[j].smiles

                row[f"reaction-product-smiles-{i}"] = reaction.product.smiles
                row[f"reaction-name-{i}"] = reaction.type
                row[f"reaction-recipe-{i}"] = None
                row[f"reaction-groupby-column-{i}"] = None
                # row[f'reaction-id-{i}'] = int(reaction.id)

            rows.append(row)

        df = DataFrame(rows)

        if len(df[df.duplicated()]):
            mrich.warning("Removing duplicates from CAR DataFrame")
            df = df.drop_duplicates()

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

        ### Get lookup data

        route_ids = self.get_routes(return_ids=True)
        smiles_lookup = self.db.get_compound_id_smiles_dict(self.reactants.compounds)
        inchikey_lookup = self.db.get_compound_id_inchikey_dict(
            self.reactants.compounds
        )

        ### Reactant Dataframe

        df = self.reactants.df

        df["smiles"] = df["compound_id"].apply(lambda x: smiles_lookup[x])
        df["inchikey"] = df["compound_id"].apply(lambda x: inchikey_lookup[x])
        df = df.drop(columns=["supplier", "max_lead_time", "quoted_amount"])

        ### Quote DataFrame

        qdf = self.db.get_quote_df(self.reactants.quote_ids)

        qdf = qdf.rename(
            columns={
                "id": "quote_id",
                "smiles": "quoted_smiles",
                "purity": "quoted_purity",
                "date": "quote_date",
                "lead_time": "quote_lead_time_days",
                "price": "quote_price",
                "currency": "quote_currency",
                "catalogue": "quote_catalogue",
                "supplier": "quote_supplier",
                "entry": "quote_entry",
                "amount": "quoted_amount_mg",
            }
        )
        qdf = qdf.drop(columns=["compound"])

        ### Join and reformat

        df = df.merge(qdf, on="quote_id", how="left")

        df = df.rename(
            columns={
                "amount": "required_amount_mg",
            }
        )

        cols = [
            "compound_id",
            "smiles",
            "inchikey",
            "required_amount_mg",
            "quoted_amount_mg",
            "quote_id",
            "quote_supplier",
            "quote_catalogue",
            "quote_entry",
            "quote_price",
            "quote_currency",
            "quote_lead_time_days",
            "quoted_purity",
            "quoted_smiles",
            "quote_date",
        ]

        df = df[[c for c in cols if c in df.columns]]

        # for i, row in mrich.track(
        #     df.iterrows(),
        #     total=len(df),
        #     prefix="Adding downstream info",
        # ):

        #     downstream_routes = []
        #     downstream_reactions = []

        #     reactant_id = row["compound_id"]

        #     for route in routes:
        #         if reactant_id in route.reactants.ids:
        #             downstream_routes.append(route)
        #         for reaction in route.reactions:
        #             if reactant_id in reaction.reactants.ids:
        #                 downstream_reactions.append(reaction)

        return df

        #     downstream_products = CompoundSet(
        #         self.db, set(route.product.id for route in downstream_routes)
        #     )
        #     downstream_reactions = ReactionSet(
        #         self.db, set(reaction.id for reaction in downstream_reactions)
        #     )

        #     if not downstream_products:
        #         mrich.error("No downstream products for", reactant)
        #         continue

        #     if not downstream_reactions:
        #         mrich.error("No downstream reactions for", reactant)
        #         continue

        #     def get_scaffold_series():

        #         bases = downstream_products.bases

        #         if not bases:
        #             bases = downstream_products[0:]

        #         return bases.ids

        #     d["num_reaction_dependencies"] = len(downstream_reactions)
        #     d["num_product_dependencies"] = len(downstream_products)
        #     d["reaction_dependencies"] = downstream_reactions.ids
        #     d["product_dependencies"] = downstream_products.ids
        #     d["chemistry_types"] = ", ".join(set(downstream_reactions.types))
        #     d["scaffold_series"] = get_scaffold_series()

        #     data.append(d)

        # df = DataFrame(data)
        # mrich.writing(file)
        # df.to_csv(file, index=False)

        # if return_df:
        #     return df

        # return None

    def write_product_csv(
        self, file: "str | Path", return_df: bool = False
    ) -> "pd.DataFrame | None":
        """Detailed CSV output including product information for selection and synthesis"""

        from pandas import DataFrame

        # from rich import print
        from .pset import PoseSet
        from .cset import CompoundSet
        from .rset import ReactionSet

        data = []

        routes = self.get_routes()

        pose_map = self.db.get_compound_id_pose_ids_dict(self.products.compounds)

        inspiration_map = self.db.get_compound_id_inspiration_ids_dict()

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

            if not upstream_routes:
                mrich.error("No upstream routes for", product)
                continue

            if not upstream_reactions:
                mrich.error("No upstream reactions for", product)
                continue

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

            inspirations = inspiration_map.get(product.id, None)

            if not inspirations and not is_base:
                base = product.bases[0]
                inspirations = inspiration_map.get(base.id, None)

                if not inspirations and "inspiration_pose_ids" in base.metadata:
                    inspirations = base.metadata["inspiration_pose_ids"]

            if (
                not inspirations
                and is_base
                and "inspiration_pose_ids" in product.metadata
            ):
                inspirations = product.metadata["inspiration_pose_ids"]

            if inspirations:
                inspirations = PoseSet(self.db, inspirations)
                d["inspirations"] = ", ".join(n for n in inspirations.names)
            else:
                d["inspirations"] = ""

            data.append(d)

        df = DataFrame(data)
        mrich.writing(file)
        df.to_csv(file, index=False)

        if return_df:
            return df

        return None

    def write_chemistry_csv(
        self, file: "str | Path", return_df: bool = True
    ) -> "pd.DataFrame | None":
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

    def to_syndirella(
        self,
        out_key: "str | Path",
        poses: "PoseSet",
        *,
        separate: bool = False,
    ) -> "DataFrame":
        """Generate inputs for running syndirella elaboration"""

        import shutil
        from pathlib import Path

        out_key = Path(".") / out_key
        out_dir = out_key.parent
        out_key = out_key.name

        mrich.var("out_key", out_key)
        mrich.var("out_dir", out_dir)

        if not out_dir.exists():
            mrich.writing(out_dir)
            out_dir.mkdir(parents=True, exist_ok=True)

        template_dir = out_dir / "templates"
        if not template_dir.exists():
            mrich.writing(template_dir)
            template_dir.mkdir(parents=True, exist_ok=True)

        """

        Need to create dataframe with columns:
        - compound_id
        - pose_id
        - smiles
        - reaction_name_step1
        - reactant_step1
        - reactant2_step1
        - product_step1
        ...
        - hit1
        - hit2
        ...
        - template
        - compound_set

        """

        pose_compounds = poses.compounds
        assert set(self.products.compound_ids) == set(
            pose_compounds.ids
        ), "supplied poses have different compounds to Recipe products"
        assert len(poses) == len(
            self.products
        ), "some duplicate compounds in supplied poses"

        df = poses.get_df(
            inchikey=False,
            alias=False,
            name=False,
            compound_id=True,
            reference_id=True,
            inspiration_aliases=True,
        )

        df = df.reset_index()
        df = df.rename(columns={"id": "pose_id"})
        df["compound_set"] = df["compound_id"].apply(lambda x: f"C{x}")
        df = df.set_index(["compound_id", "pose_id"])

        ## CHECKS

        no_refs = df[df["reference_id"].isna()]

        if len(no_refs):
            mrich.error(len(no_refs), "poses without reference!")
            ids = set(no_refs.index.get_level_values("pose_id"))
            mrich.print(ids)
            from .pset import PoseSet

        no_insps = bool([1 for i in df["inspiration_aliases"].values if not len(i)])

        if no_insps:
            mrich.error(len(no_insps), "poses without inspirations!")
            return None

        ## TEMPLATES

        references = poses.references
        ref_lookup = self.db.get_pose_id_alias_dict(references)
        df["template"] = df["reference_id"].apply(lambda x: ref_lookup[x])

        for ref_pose in references:
            assert ref_pose.apo_path, f"Reference {ref_pose} has no apo_path"

            template = template_dir / ref_pose.apo_path.name

            if not template.exists():
                mrich.writing(template)
                shutil.copy(ref_pose.apo_path, template)

        ## INSPIRATIONS

        for i, row in df.iterrows():
            for j, alias in enumerate(row["inspiration_aliases"]):
                df.loc[i, f"hit{j+1}"] = alias

        inspirations = poses.inspirations

        sdf_name = out_dir / f"{out_key}_syndirella_inspiration_hits.sdf"

        inspirations.write_sdf(
            sdf_name,
            tags=False,
            metadata=False,
            name_col="name",
        )

        ## ADD ROUTE INFO

        routes = self.get_routes()

        for sub_recipe in mrich.track(routes, prefix="Adding chemistry info..."):

            product = sub_recipe.product

            product_id = product.compound_id

            matches = df.xs(product_id, level="compound_id")

            if len(matches) > 1:
                mrich.warning("Multiple rows for compound", product_id)

            for i, row in matches.iterrows():

                key = (product_id, i)

                for j, reaction in enumerate(sub_recipe.reactions):

                    j = j + 1

                    match len(reaction.reactants):
                        case 1:
                            df.loc[key, f"reactant_step{j}"] = reaction.reactants[
                                0
                            ].smiles
                            df.loc[key, f"reactant2_step{j}"] = None
                        case 2:
                            df.loc[key, f"reactant_step{j}"] = reaction.reactants[
                                0
                            ].smiles
                            df.loc[key, f"reactant2_step{j}"] = reaction.reactants[
                                1
                            ].smiles
                        case 3:
                            df.loc[key, f"reactant_step{j}"] = reaction.reactants[
                                0
                            ].smiles
                            df.loc[key, f"reactant2_step{j}"] = reaction.reactants[
                                1
                            ].smiles
                            df.loc[key, f"reactant3_step{j}"] = reaction.reactants[
                                2
                            ].smiles
                        case _:
                            raise NotImplementedError("Too many reactants")

                    df.loc[key, f"product_step1{j}"] = reaction.product.smiles
                    df.loc[key, f"reaction_name_step{j}"] = reaction.type

                break

        ## REMOVE UNECESSARY COLS

        df = df.drop(columns=["reference_id", "inspiration_aliases"])

        ## REORDER COLUMNS

        cols = [
            "smiles",
            "reaction_name_step1",
            "reactant_step1",
            "reactant2_step1",
            "reactant3_step1",
            "product_step11",
            "hit1",
            "hit2",
            "hit3",
            "hit4",
            "hit5",
            "hit6",
            "hit7",
            "hit8",
            "hit9",
            "template",
            "compound_set",
        ]

        if not any([c not in cols for c in df.columns]):
            df = df[[c for c in cols if c in df.columns]]

        if not separate:
            out_path = out_dir / f"{out_key}_syndirella_input.csv"
            mrich.writing(out_path)
            df.to_csv(out_path)
            return df

        for idx, row in df.iterrows():
            out_path = out_dir / f"{out_key}_{row['compound_set']}_syndirella_input.csv"
            mrich.writing(out_path)
            single_df = row.to_frame().T
            single_df = single_df.dropna(axis=1, how="all")
            single_df.to_csv(out_path, index=False)

        return df

    def copy(self) -> "Recipe":
        """Copy this recipe"""

        if hasattr(self, "compounds"):
            compounds = self.compounds.copy()
        else:
            compounds = None

        return Recipe(
            self.db,
            products=self.products.copy(),
            reactants=self.reactants.copy(),
            intermediates=self.intermediates.copy(),
            reactions=self.reactions.copy(),
            compounds=compounds,
            # supplier=self.supplier
        )

    def __flag_modification(self) -> None:
        """Flag this recipe as modified"""
        self._product_interactions = None
        self._score = None
        self._product_compounds = None
        self._product_poses = None

    def check_integrity(self, debug: bool = False) -> bool:
        """Verify integrity of this recipe"""

        # no duplicate ingredients

        if debug:
            mrich.debug("Checking integrity:", self)
            mrich.debug("Checking for duplicate compounds")

        if len(self.reactants.compound_ids) != len(set(self.reactants.compound_ids)):
            mrich.error("Reactant compound ID's are not unique")
            return False
        if len(self.intermediates.compound_ids) != len(
            set(self.intermediates.compound_ids)
        ):
            mrich.error("Intermediate compound ID's are not unique")
            return False
        if len(self.products.compound_ids) != len(set(self.products.compound_ids)):
            mrich.error("Product compound ID's are not unique")
            return False

        # all references should exist

        if debug:
            mrich.debug("Checking for missing references")

        if self.db.count_where(
            table="reaction", key=f"reaction_id IN {self.reactions.str_ids}"
        ) < len(self.reactions):
            mrich.error("Not all Reactions in Database")
            return False

        if self.db.count_where(
            table="compound", key=f"compound_id IN {self.product_compounds.str_ids}"
        ) < len(self.products):
            mrich.error("Not all product Compounds in Database")
            return False

        if self.db.count_where(
            table="compound", key=f"compound_id IN {self.reactants.compounds.str_ids}"
        ) < len(self.reactants):
            mrich.error("Not all reactant Compounds in Database")
            return False

        if self.db.count_where(
            table="compound",
            key=f"compound_id IN {self.intermediates.compounds.str_ids}",
        ) < len(self.intermediates):
            mrich.error("Not all intermediate Compounds in Database")
            return False

        reaction_intermediates = self.reactions.intermediates
        reaction_products = self.reactions.products
        reaction_reactants = self.reactions.reactants

        if debug:
            mrich.debug("Checking for missing reactions")

        # all products should have a reaction
        for product in self.products:
            if product not in reaction_products:
                mrich.error(f"Product: {product} does not have associated reaction")
                return False

        # intermediates
        for intermediate in self.intermediates:
            if intermediate not in reaction_intermediates:
                mrich.error(
                    f"Intermediate: {intermediate} is not in self.reactions.intermediates"
                )
                return False

        # reactants
        for reactant in self.reactants:
            if reactant not in reaction_reactants:
                mrich.error(f"Reactant: {reactant} is not in self.reactions.reactants")
                return False

        # all reactions should have enough reactant

        if debug:
            mrich.debug("Checking reactant quantities")

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

        if debug:
            mrich.success(self, "OK")

        return True

    def add_ingredient(self, ingredient: "Ingredient", amount: float = 1):
        """Add an :class:`.Ingredient` object for direct purchase (no associated reactions)"""
        self.compounds.add(ingredient)

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""

        if self.score:
            s = f"(score={self.score:.3f})"
        else:
            s = ""

        if self.hash:
            return f"Recipe_{self.hash}{s}"

        return f"Recipe{s}"

    def __longstr(self) -> str:
        """Unformatted string representation"""

        if self.empty:
            return f"Empty Recipe()"

        if self.reactions:

            if self.intermediates:
                s = f"{self.reactants} --> {self.intermediates} --> {self.products} via {self.reactions}"
            else:
                s = f"{self.reactants} --> {self.products} via {self.reactions}"

            if self.score:
                s += f", score={self.score:.3f}"

            if self.hash:
                return f"Recipe_{self.hash}({s})"

            return f"Recipe({s})"

        else:

            s = f"{self.compounds}"

            if self.hash:
                return f"Recipe_{self.hash}({s})"

            return f"Recipe(#compounds={self.num_compounds} [no-chem])"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        return f"{mcol.bold}{mcol.underline}{self.__longstr()}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self.__longstr()}"

    def __add__(self, other: "Recipe"):
        result = self.copy()
        result.reactants += other.reactants
        result.intermediates += other.intermediates
        result.reactions += other.reactions
        result.products += other.products
        if hasattr(other, "compounds"):
            result.compounds += other.compounds
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
    def from_json(
        cls, db: "Database", path: "str | Path", data: dict = None
    ) -> "Route":
        """Load a serialised route from a JSON file

        :param db: database to link
        :param path: path to JSON
        :param data: serialised data (Default value = None)

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
    def product(self) -> "Ingredient":
        """Product ingredient"""
        return self._products[0]

    @property
    def product_compound(self) -> "Compound":
        """Product compound"""
        return self.product.compound

    @property
    def id(self) -> int:
        """Route ID"""
        return self._id

    @property
    def price(self) -> "Price":
        """Get the price of the reactants"""
        return self.reactants.price

    ### METHODS

    def get_dict(self) -> dict:
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
    def from_json(
        cls, db: "Database", path: "str | Path", data: dict = None
    ) -> "RouteSet":
        """Load a serialised routeset from a JSON file

        :param db: database to link
        :param path: path to JSON
        :param data: serialised data (Default value = None)

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
    def data(self) -> "dict[int, Route]":
        """Get internal data dictionary"""
        return self._data

    @property
    def db(self):
        """Get associated database"""
        return self._db

    @property
    def routes(self) -> "list[Route]":
        """Get route objects"""
        return self.data.values()

    @property
    def product_ids(self) -> list[int]:
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
    def str_ids(self) -> str:
        """Return an SQL formatted tuple string of the :class:`.Route` ID's"""
        return str(tuple(self.ids)).replace(",)", ")")

    @property
    def ids(self) -> list[int]:
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

    def copy(self) -> "RouteSet":
        """Copy this RouteSet"""
        return RouteSet(self.db, self.data.values())

    def set_db_pointers(self, db: "Database") -> None:
        """

        :param db:

        """
        self._db = db
        for route in self.data.values():
            route._db = db

    # def clear_db_pointers(self):
    #     """ """
    #     self._db = None
    #     for route in self.data.values():
    #         route._db = None

    # def get_dict(self):
    #     """Get serialisable dictionary"""

    #     data = dict(db=str(self.db), routes={})

    #     # populate with routes
    #     for route_id, route in self.data.items():
    #         data["routes"][route_id] = route.get_dict()

    #     return data

    def pop_id(self) -> int:
        """Pop the last route from the set and return it's id"""
        route_id, route = self.data.popitem()
        return route_id

    def pop(self) -> "Route":
        """Pop the last route from the set and return it's object"""
        route_id, route = self.data.popitem()
        return route

    def balanced_pop(
        self, permitted_clusters: set[tuple] | None = None, debug: bool = False
    ) -> "Route":
        """Pop a route from this set, while maintaining the balance of scaffold clusters populations"""

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

    def __len__(self) -> int:
        """Number of routes in this set"""
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
    """A set of recipes stored on disk"""

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
    def db(self) -> "Database":
        """Associated database"""
        return self._db

    ### METHODS

    def get_values(
        self,
        key: str,
        progress: bool = False,
        serialise_price: bool = False,
    ):
        """Get values of member recipes associated with attribute ``key``

        :param key: attribute to query/calculate
        :param progress: show a progress bar
        :param serialise_price: serialise price objects

        """

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
        """Get dataframe of recipe dictionaries. See :meth:`.Recipe.get_dict`"""

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

    def items(self) -> "list[tuple[str, Recipe]]":
        """Get data dictionary items"""
        return self._recipes.items()

    def keys(self) -> list[str]:
        """Get data dictionary keys (recipe hashes)"""
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
