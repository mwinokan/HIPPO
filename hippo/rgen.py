from .recipe import Recipe
from .cset import CompoundSet, IngredientSet

from tqdm import tqdm

from pathlib import Path
import json

from .tools import dt_hash

import mrich


class RandomRecipeGenerator:
    """ """

    """Class to create randomly sampled Recipe from a HIPPO Database"""

    def __init__(
        self,
        db,
        *,
        max_lead_time=None,
        # max_reactions = None,
        suppliers: list | None = None,
        start_with: Recipe | CompoundSet | IngredientSet = None,
    ):

        mrich.debug("RandomRecipeGenerator.__init__()")

        # Static parameters
        self._db_path = db.path
        self._max_lead_time = max_lead_time
        self._suppliers = suppliers
        self._starting_recipe = start_with

        mrich.var("database", self.db_path)
        mrich.var("max_lead_time", self.max_lead_time)
        mrich.var("suppliers", self.suppliers)

        # Database set up
        self._db = db

        # JSON I/O set up
        self._data_path = Path(str(self.db_path).replace(".sqlite", "_rgen.json"))
        if self.data_path.exists():
            mrich.warning(f"Will overwrite existing rgen data file: {self.data_path}")

        # Recipe I/O set up
        path = Path(str(self.db_path).replace(".sqlite", "_recipes"))
        mrich.writing(f"{path}/")
        path.mkdir(exist_ok=True)
        self._recipe_dir = path

        # Route pool
        mrich.debug("Solving route pool...")
        self._route_pool = self.get_route_pool()

        # dump data
        self.dump_data()

    ### FACTORIES

    @classmethod
    def from_json(cls, db: "Database", path: "Path | str"):
        """Construct the RandomRecipeGenerator from a JSON file

        :param db:
        :param path:

        """

        data = json.load(open(path, "rt"))

        self = cls.__new__(cls)

        self._db_path = Path(data["db_path"])
        self._recipe_dir = Path(data["recipe_dir"])
        self._max_lead_time = data["max_lead_time"]
        self._suppliers = data["suppliers"]

        self._starting_recipe = Recipe.from_json(
            db=db,
            path=None,
            data=data["starting_recipe"],
            allow_db_mismatch=True,
        )

        mrich.var("database", self.db_path)
        mrich.var("max_lead_time", self.max_lead_time)
        mrich.var("suppliers", self.suppliers)

        self._db = db

        # JSON I/O set up
        self._data_path = Path(path)

        # Route pool
        from .recipe import RouteSet

        self._route_pool = RouteSet.from_json(path=None, data=data["route_pool"], db=db)

        return self

    ### PROPERTIES

    @property
    def starting_recipe(self):
        """Get the starting recipe used in all generations"""
        return self._starting_recipe

    @property
    def db(self) -> "Database":
        """Get the linked HIPPO Database object"""
        return self._db

    @property
    def db_path(self) -> str:
        """Get the path of the linked Database"""
        return self._db_path

    @property
    def suppliers_str(self) -> str:
        """SQL formatted tuple of suppliers"""
        return str(tuple(self.suppliers)).replace(",)", ")")

    @property
    def suppliers(self) -> list[str]:
        """List of suppliers"""
        return self._suppliers

    @property
    def max_lead_time(self) -> float:
        """Maximum lead-time constraint"""
        return self._max_lead_time

    @property
    def route_pool(self):
        """Get the RouteSet of all product reaction routes considered by this generator"""
        return self._route_pool

    @property
    def data_path(self):
        """File path for the JSON data export"""
        return self._data_path

    @property
    def recipe_dir(self):
        """File path for the JSON recipe export"""
        return self._recipe_dir

    ### POOL METHODS

    def get_route_pool(self, mini_test=False):
        """Construct the pool of routes that will be randomly sampled from

        :param mini_test:  (Default value = False)

        """

        """
			Explainer for SQL query:

			- get table of quoted compounds with a count of the valid suppliers
			- join routes, components, and the new table together and grouped by route count the unavailable reactants
			- return route ids where no reactants are unavailable

		"""

        if "route" not in self.db.table_names:
            mrich.error("route table not in Database")
            raise NotImplementedError

        assert self.suppliers_str
        if self.max_lead_time:
            raise NotImplementedError

        ### EXCLUDE PRODUCTS OF ROUTES IN STARTING RECIPE!!!

        route_ids = self.db.execute(
            f"""
		WITH possible_reactants AS (
			SELECT quote_compound, COUNT(CASE WHEN quote_supplier IN {self.suppliers_str} THEN 1 END) AS [count_valid] FROM quote
			GROUP BY quote_compound
		),

		route_reactants AS (
			SELECT route_id, route_product, 
            COUNT(
                CASE 
                    WHEN count_valid = 0 THEN 1 
                    WHEN count_valid IS NULL THEN 1 
                END) 
            AS [count_unavailable] FROM route
			INNER JOIN component ON component_route = route_id
			LEFT JOIN possible_reactants ON quote_compound = component_ref
			WHERE component_type = 2
			GROUP BY route_id
		)

		SELECT route_id FROM route_reactants
		WHERE count_unavailable = 0
		"""
        ).fetchall()

        if mini_test:
            route_ids = route_ids[:100]

        routes = [
            self.db.get_route(id=route_id)
            for route_id, in mrich.track(route_ids, prefix="Getting routes")
        ]

        from .recipe import RouteSet

        return RouteSet(self.db, routes)

    ### FILE I/O METHODS

    def dump_data(self):
        """ """

        data = {}

        data["db_path"] = str(self.db_path.resolve())
        data["recipe_dir"] = str(self.recipe_dir.resolve())
        data["max_lead_time"] = self.max_lead_time
        data["suppliers"] = self.suppliers
        data["starting_recipe"] = self.starting_recipe.get_dict(serialise_price=True)
        data["route_pool"] = self.route_pool.get_dict()

        mrich.writing(self.data_path)
        json.dump(data, open(self.data_path, "wt"), indent=4)

    def generate(
        self,
        budget: float = 10000,
        currency: str = "EUR",
        max_products=1000,
        max_reactions=1000,
        debug=False,
        max_iter=None,
        shuffle=True,
        balance_clusters: bool = False,
        permitted_clusters: None | set = None,
    ):
        """

        :param budget: float:  (Default value = 10000)
        :param currency: str:  (Default value = 'EUR')
        :param max_products:  (Default value = 1000)
        :param max_reactions:  (Default value = 1000)
        :param debug:  (Default value = True)
        :param max_iter:  (Default value = None)
        :param # pick_inner_cheapest:  (Default value = True)
        :param # add_size:  (Default value = 1)
        :param shuffle:  (Default value = True)

        """

        # construct filename

        out_file = self.recipe_dir / f"Recipe_{dt_hash()}.json"

        from .price import Price

        if not max_iter:
            max_iter = max_products + max_reactions

        budget = Price(budget, currency)

        recipe = self.starting_recipe.copy()

        recipe.reactants._supplier = self.suppliers

        # get the RouteSet
        pool = self.route_pool.copy()

        assert len(pool), "Route pool is empty!"

        if shuffle:
            mrich.debug("Shuffling Route pool")
            pool.shuffle()

        old_recipe = recipe.copy()

        mrich.var("route pool", len(pool))
        mrich.var("max_iter", max_iter)

        for i in mrich.track(range(max_iter), prefix="Generating Recipe..."):

            if debug:
                mrich.title(f"Iteration {i}")

            price = recipe.price
            mrich.set_progress_field("price", str(price))
            mrich.set_progress_field("#products", len(recipe.products))

            if debug:
                mrich.var("price", price)

            # pop a route
            if balance_clusters:
                candidate_route = pool.balanced_pop(
                    permitted_clusters=permitted_clusters
                )
            else:
                candidate_route = pool.pop()

            if debug:
                mrich.var("candidate_route", candidate_route)
            if debug:
                mrich.var("candidate_route.reactants", candidate_route.reactants.ids)

            if candidate_route.product in recipe.products:
                continue

            # add the route to the recipe
            if debug:
                mrich.var("#recipe.reactants", len(recipe.reactants))
            recipe += candidate_route
            if debug:
                mrich.var("#recipe.reactants", len(recipe.reactants))

            # calculate the new price
            try:
                new_price = recipe.price
            except AssertionError:
                mrich.error(
                    f"Something went wrong while calculating the price after adding {candidate_route=} to recipe"
                )
                raise

            if debug:
                mrich.var("new price", new_price)

            # Break if product pool depleted
            if not len(pool):
                stop_reason = "Product pool depleted"
                mrich.success(stop_reason)
                break

            # check breaking conditions
            if new_price > budget:
                recipe = old_recipe.copy()
                continue

            if len(recipe.reactions) > max_reactions:
                stop_reason = "Max #reactions exceeded"
                mrich.success(stop_reason)
                break

            if len(recipe.products) > max_products:
                stop_reason = "Max #products exceeded"
                mrich.success(stop_reason)
                break

            # accept change
            old_recipe = recipe.copy()

        else:
            stop_reason = "Max #iterations reached"
            mrich.warning(stop_reason)

        ### recalculate the products to see if any extra can be had for free?

        mrich.success(f"Completed after {i} iterations")

        metadict = {
            "rgen_data_path": str(self.data_path.resolve()),
            "rgen_db_path": str(self.db_path.resolve()),
            "rgen_recipe_dir": str(self.recipe_dir.resolve()),
            "rgen_max_lead_time": self.max_lead_time,
            "rgen_suppliers": self.suppliers,
            "gen_budget": budget.amount,
            "gen_currency": budget.currency,
            "gen_max_products": max_products,
            "gen_max_reactions": max_reactions,
            "gen_max_iter": max_iter,
            "gen_shuffle": shuffle,
            "gen_iterations": i,
            "gen_stop_reason": stop_reason,
            "gen_recipe_path": str(out_file.resolve()),
        }

        # write the Recipe JSON
        recipe.write_json(out_file, extra=metadict)

        return recipe

    ### DUNDERS

    def __call__(self, *args, **kwargs):
        return self.generate(*args, **kwargs)

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"RandomRecipeGenerator(recipe_dir={self.recipe_dir})"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        import mcol

        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"
