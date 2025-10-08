import mout
import numpy as np

# from .block import BuildingBlockSet
from scipy.interpolate import interp1d
import pandas as pd

import mrich

DATA_COLUMNS = [
    "score",
    "price",
    "compound_ids",
    "pose_ids",
    "interaction_ids",
]


class Scorer:
    """

    :param recipes: :class:`.RecipeSet`

    inputs
    ------

    * building block set
    * HIPPO object

    parameters
    ----------

    * optimisation attributes
    * weights

    outputs
    -------

    * score for a given BBS
    * BBSs sorted by score

    """

    def __init__(
        self,
        db: "Database",
        directory: "Path | str",
        pattern: str = "*.json",
        attributes: list[str] = None,
        populate: bool = True,
        load_cache: bool = True,
    ) -> None:

        from .recipe import RecipeSet

        self._db = db

        attributes = attributes or []

        rset = RecipeSet(db, directory, pattern=pattern)

        self._recipes = rset

        self._attributes = {}

        for key in attributes:
            attribute = Attribute(self, key)
            self._attributes[key] = attribute

        self._data = pd.DataFrame(
            index=recipes.keys(),
            columns=DATA_COLUMNS + self.attribute_keys,
        )

        self._data.replace({np.nan: None}, inplace=True),

        if populate:
            if load_cache and self.json_path.exists():
                self._load_json()
            else:
                self._populate_query_cache()

            self._populate_recipe_child_sets()

        self.weights = 1.0

    ### FACTORIES

    @classmethod
    def default(
        cls,
        db: "Database",
        directory: "Path | str",
        pattern: str = "*.json",
        skip: list[str] | None = None,
        load_cache: bool = True,
        subsites: bool = True,
    ) -> "Scorer":
        """Create a Scorer instance with Default attributes"""

        from .recipe import RecipeSet

        self = cls.__new__(cls)

        attributes = [
            k for k, v in DEFAULT_ATTRIBUTES.items() if v["type"] == "standard"
        ]

        self.__init__(
            db=db,
            directory=directory,
            pattern=pattern,
            attributes=attributes,
            populate=False,
        )

        skip = skip or []

        if not db.count("interaction"):
            mrich.warning("No interactions in DB, skipping related metrics")
            skip.append("interaction_count")
            skip.append("interaction_balance")

        if not db.count("pose"):
            mrich.warning("No poses in DB, skipping related metrics")
            skip.append("num_inspirations")
            skip.append("num_inspiration_sets")
            skip.append("avg_energy_score")
            skip.append("avg_distance_score")

        if not db.count("scaffold"):
            mrich.warning("No scaffold entries in DB, skipping related metrics")
            skip.append("num_scaffolds")
            skip.append("num_scaffolds_elaborated")
            skip.append("elaboration_balance")

        # custom attributes
        for key, attribute in [
            (k, v) for k, v in DEFAULT_ATTRIBUTES.items() if v["type"] == "custom"
        ]:

            if skip and key in skip:
                continue

            if not subsites and "subsite" in key:
                continue

            self.add_custom_attribute(
                key, attribute["function"], weight_reset_warning=False
            )

        if load_cache and self.json_path.exists():
            self._load_json()
        else:
            self._populate_query_cache()

        self._populate_recipe_child_sets()

        # weights
        wsum = sum(abs(d["weight"]) for d in DEFAULT_ATTRIBUTES.values())
        for attribute in self.attributes:
            d = DEFAULT_ATTRIBUTES[attribute.key]
            attribute.weight = d["weight"] / wsum

        return self

    ### PROPERTIES

    @property
    def num_recipes(self) -> int:
        """Number of recipes being evaluated"""
        return len(self._recipes)

    @property
    def attributes(self):
        return list(self._attributes.values())

    @property
    def attribute_keys(self):
        return list(self._attributes.keys())

    @property
    def recipes(self):
        return self._recipes

    @property
    def num_attributes(self):
        return len(self.attributes)

    @property
    def weights(self):
        return [a.weight for a in self.attributes]

    @weights.setter
    def weights(self, ws):

        self._flag_weight_modification()

        if isinstance(ws, float):
            ws = [ws] * self.num_attributes

        ws = [w for w in ws]
        wsum = sum([abs(w) for w in ws])

        for a, w in zip(self.attributes, ws):
            a.weight = w / wsum

    @property
    def score_dict(self):

        col = self._data["score"]

        null = col.isnull()

        if null.sum():
            mrich.debug("Calculating scores...")
            for key in col[null].index.values:
                recipe = self.recipes[key]
                score = self.score(recipe)
                self._data.at[key, "score"] = score

            self._dump_json()

        return col.to_dict()

    @property
    def scores(self):
        return list(self.score_dict.values())

    @property
    def best(self):
        return self.top(1)

    @property
    def db(self):
        return self._db

    @property
    def json_path(self):
        from pathlib import Path

        return Path(self.db.path.name.replace(".sqlite", "_scorer.json"))

    ### METHODS

    def add_custom_attribute(
        self,
        key,
        function,
        weight_reset_warning: bool = True,
    ) -> "CustomAttribute":

        ca = CustomAttribute(self, key, function)

        if key not in self._attributes:

            # self._flag_weight_modification()

            self._attributes[key] = ca

            if weight_reset_warning:
                mrich.warning("Attribute weights have been reset")
            self.weights = 1.0

            self._data[key] = None

        else:
            mrich.warning("Existing attribute with {key=}")

        return self._attributes[key]

    def add_recipes(self, json_paths: "list", debug: bool = False):

        from pathlib import Path
        from .recipe import Recipe

        for json_path in json_paths:

            path = Path(json_path)

            key = path.name.removeprefix("Recipe_").removesuffix(".json")

            if key in self.recipes:
                mrich.warning(f"Skipping duplicate {path}")
                continue

            recipe = Recipe.from_json(self._db, path, allow_db_mismatch=True)

            recipe._hash = key

            if debug:
                mrich.debug(recipe)

            if debug:
                mrich.debug("Updating Scorer.recipes._json_paths")
            self.recipes._json_paths[key] = path.resolve()

            if debug:
                mrich.debug("Updating Scorer.recipes._recipes")
            self.recipes._recipes[key] = recipe

            self._data.loc[key] = None

        self._data.replace({np.nan: None}, inplace=True),
        self._populate_query_cache()
        self._populate_recipe_child_sets()
        self._flag_weight_modification()

    def score(
        self,
        recipe: "Recipe",
        *,
        debug: bool = False,
    ) -> float:

        score = 0.0

        for attribute in self.attributes:
            recipe_score = attribute(recipe)
            score += recipe_score
            # score = sum([a(recipe) for a in self.attributes])

        if debug:
            print_data = []
            for attribute in self.attributes:
                print_data.append(
                    dict(
                        key=attribute.key,
                        weight=attribute.weight,
                        value=attribute.get_value(recipe),
                        unweighted=attribute.unweighted(recipe),
                        weighted=attribute(recipe),
                    )
                )

            print(pd.DataFrame(print_data))
            mrich.var("score", score)

        recipe._score = score

        return score

    def get_df(
        self,
        serialise_price: bool = True,
        debug: bool = True,
        **kwargs,
    ) -> "pandas.DataFrame":

        raise NotImplementedError

        from pandas import DataFrame

        params = dict(serialise_price=serialise_price, **kwargs)

        if self._df is None or params != self._df_params:

            if debug:
                mrich.debug("Scorer.get_df()")

            if debug:
                mrich.debug("Scorer.recipes.get_df()")

            data = []

            for recipe in self.recipes:

                key = recipe.hash

                d = recipe.get_dict(
                    # reactant_supplier=False,
                    database=False,
                    timestamp=False,
                    **kwargs,
                    # timestamp=False,
                )

                d["hash"] = key
                d["score"] = self.scores[key]

                data.append(d)

            df = DataFrame(data)

            for attribute in self.attributes:

                if debug:
                    mrich.debug(f"Getting {attribute.key=} values")

                if isinstance(attribute, CustomAttribute):

                    df[attribute.key] = attribute.values

                else:

                    df[attribute.key] = self.recipes.get_values(
                        key=attribute.key, serialise_price=serialise_price
                    )

            if serialise_price:
                df["price"] = df.apply(lambda x: x["price"].amount, axis=1)

            self._df = df
            self._df_params = params

        return self._df

    def get_sorted_df(self, budget: float | None = None):
        self.scores
        return self._data.sort_values(by="score", ascending=False)

    def plot(
        self,
        keys: str | list[str],
        budget: float | None = None,
        # **kwargs,
    ):

        import plotly.express as px

        if len(keys) != 2:
            mrich.error("Only two keys supported")
            return None

        # calculate scores
        self.scores

        df = self._data.drop(
            columns=[
                "compound_ids",
                "pose_ids",
                "interaction_ids",
            ]
        )

        df["score"] = pd.to_numeric(df["score"])

        if isinstance(keys, str):
            assert keys in df.columns
            return px.histogram(df, x=keys)

        if not all(key in df.columns for key in keys):
            for key in keys:
                if key not in df.columns:
                    raise KeyError(f'no attribute/column named "{key}"')

        if budget:
            df = df[df["price"] < budget]

        df["hash"] = df.index.values

        hover_data = [
            "hash",
        ]

        hover_data += [c for c in df.columns]

        return px.scatter(
            df, x=keys[0], y=keys[1], color="score", hover_data=hover_data
        )

    def top_keys(self, n: int, budget: float | None = None):
        keys = self.get_sorted_df(budget=budget)["hash"][:n]
        return keys

    def top(self, n: int, budget: float | None = None):
        keys = self.top_keys(n=n, budget=budget)
        if n == 1:
            return [self.recipes[key] for key in keys][0]
        else:
            return [self.recipes[key] for key in keys]

    ### INTERNALS

    def _flag_weight_modification(self):

        self._data["score"] = None

    def summary(self):

        mrich.header(self)
        for attribute in self.attributes:
            # mrich.out(
            # f"{repr(attribute)} min={attribute.min:.3g}, mean={attribute.mean:.3g}, std={attribute.std:.3g}, max={attribute.max:.3g}"
            # )
            mrich.print(
                attribute,
                f"min={attribute.min:.3g}, mean={attribute.mean:.3g}, std={attribute.std:.3g}, max={attribute.max:.3g}",
            )

    def __check_integrity(self):

        n_recipes = len(self.recipes)

        for attribute in self.attributes:
            assert len(attribute._value_dict) == n_recipes

        assert len(self._scores) == n_recipes

        assert len(self._data) == n_recipes
        assert len(self._data.columns) == len(attributes) + len(DATA_COLUMNS)

        return True

    def _populate_query_cache(self):

        from .cset import CompoundSet
        from .pset import PoseSet

        df = self._data

        ### Recipe prices

        for recipe in self.recipes:
            self._data.at[recipe.hash, "price"] = recipe.price.amount

        ### Product Compound IDs

        col = "compound_ids"
        null = df[col].isnull()

        # populate missing product compound ids
        if null.sum():
            mrich.debug(f'Populating _data["{col}"]...')
            assert len(df[null]) == null.sum()
            for key in df[null].index.values:
                recipe = self.recipes[key]
                df.at[key, col] = recipe.combined_compound_ids

        ### Product Pose IDs

        col = "pose_ids"
        null = df[col].isnull()

        # populate missing product pose ids
        if null.sum():

            compound_ids = set()
            for ids in df[null]["compound_ids"]:
                for id in ids:
                    compound_ids.add(id)

            cset = CompoundSet(self.db, compound_ids, sort=False)

            mrich.debug(f"Getting poses for {len(cset)} compounds")
            pose_map = self.db.get_compound_id_pose_ids_dict(cset)

            mrich.debug(f'Populating _data["{col}"]...')
            for key in df[null].index.values:
                assert len(df[null]) == null.sum()
                recipe = self.recipes[key]
                comp_ids = df["compound_ids"][key]

                row = df.loc[key]

                all_pose_ids = set()

                for comp_id in comp_ids:
                    pose_ids = pose_map.get(comp_id, set())
                    all_pose_ids |= pose_ids

                df.at[key, col] = all_pose_ids

        ### Product Interaction IDs

        col = "interaction_ids"
        null = df[col].isnull()

        # populate missing product interaction ids
        if null.sum():

            pose_ids = set()
            for ids in df[null]["pose_ids"]:
                for id in ids:
                    pose_ids.add(id)

            pset = PoseSet(self.db, pose_ids, sort=False)

            mrich.debug(f"Getting interactions for {len(pset)} poses")
            interaction_map = self.db.get_pose_id_interaction_ids_dict(pset)

            mrich.debug(f'Populating _data["{col}"]...')
            for key in df[null].index.values:
                assert len(df[null]) == null.sum()
                recipe = self.recipes[key]
                pose_ids = df["pose_ids"][key]

                row = df.loc[key]

                all_interaction_ids = set()

                for pose_id in pose_ids:
                    interaction_ids = interaction_map.get(pose_id, set())
                    all_interaction_ids |= interaction_ids

                df.at[key, col] = all_interaction_ids

        # raise NotImplementedError

    def _populate_recipe_child_sets(self):

        from .cset import CompoundSet
        from .pset import PoseSet
        from .iset import InteractionSet

        mrich.debug("Populating recipe caches")
        for key, recipe in self.recipes.items():

            row = self._data.loc[key]

            if recipe._combined_compounds is None:
                ids = row["compound_ids"]
                cache = CompoundSet(self.db, ids)
                cache._name = f"Recipe_{key} products"
                recipe._combined_compounds = cache

            if recipe._poses is None:
                ids = row["pose_ids"]
                cache = PoseSet(self.db, ids)
                cache._name = f"Recipe_{key} product poses"
                recipe._product_poses = cache

            if recipe._interactions is None:
                ids = row["interaction_ids"]
                cache = InteractionSet(self.db, ids)
                cache._name = f"Recipe_{key} product interactions"
                recipe._interactions = cache

    def _dump_json(self):
        path = self.json_path
        mrich.writing(path)
        self._data.to_json(path)

    def _load_json(self):
        path = self.json_path

        mrich.reading(path)
        cached = pd.read_json(path, orient="columns")

        if (cached_columns := set(cached.columns)) != (
            self_columns := set(self._data.columns)
        ):

            for col in cached_columns - self_columns:
                mrich.error(f"JSON has unexpected {col}")

            for col in self_columns - cached_columns:
                mrich.error(f"JSON is missing {col}")

            display(cached.head())
            display(self._data.head())

            raise ValueError("JSON columns don't match expectation")

        cached_keys = set(cached.index.values)
        self_keys = set(self._data.index.values)

        if difference := cached_keys - self_keys:
            mrich.warning("JSON has extra Recipes:")
            mrich.warning(difference)

        if difference := self_keys - cached_keys:
            mrich.error("JSON is missing Recipes:")
            mrich.error(difference)
            raise ValueError("JSON is missing Recipes")

        cached.replace({np.nan: None}, inplace=True),

        self._data = cached

        # return cached

    ### DUNDERS

    def __str__(self) -> str:
        """Unformatted string representation"""
        return f"Scorer(#recipes={self.num_recipes})"

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        import mcol

        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"


class Attribute:

    _type = "Attribute"

    ### DUNDERS

    def __init__(
        self,
        scorer: "Scorer",
        key: str,
        *,
        inverse: bool = False,
        weight: float = 1.0,
        bins: int = 100,
    ):

        self._scorer = scorer

        self._key = key
        self._inverse = inverse
        self._weight = weight

        self._value_dict = {}

        self._bins = bins

        self._percentile_interpolator = None

    ### PROPERTIES

    @property
    def scorer(self):
        return self._scorer

    @property
    def key(self):
        return self._key

    @property
    def inverse(self):
        return self._inverse

    @property
    def bins(self):
        return self._bins

    @property
    def value_dict(self):
        df = self.scorer._data[self.key]

        null = df.isnull()

        if null.sum():
            # for key in mrich.track(
            # df[null].index.values, f"Constructing value dictionary for {self}"
            # , total = len(df[null])):
            with mrich.loading(f"Constructing value dictionary for {self}"):
                for key in df[null].index.values:
                    recipe = self.scorer.recipes[key]
                    self.get_value(recipe, force=True)
                self.scorer._dump_json()

        return df.to_dict()

    @property
    def values(self):
        return list(self.value_dict.values())

    @property
    def mean(self):
        return np.mean(self.values)

    @property
    def std(self):
        return np.std(self.values)

    @property
    def max(self):
        return max(self.values)

    @property
    def min(self):
        return min(self.values)

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, w):
        self.scorer._flag_weight_modification()
        self._weight = abs(w)
        self._reverse = w < 0

    @property
    def percentile_interpolator(self):
        if self._percentile_interpolator is None:

            count, bins_count = np.histogram(self.values, bins=self.bins)

            pdf = count / sum(count)
            cdf = np.cumsum(pdf)
            self._percentile_interpolator = interp1d(
                bins_count[1:], cdf, kind="linear", fill_value="extrapolate"
            )

        return self._percentile_interpolator

    ### METHODS

    def get_value(
        self,
        recipe: "Recipe",
        serialise_price: bool = True,
        force: bool = False,
    ) -> float:

        if not force:
            cached = self.scorer._data[self.key][recipe.hash]

        if force or cached is None:
            value = getattr(recipe, self.key)
            if serialise_price and self.key == "price":
                value = value.amount
            self.scorer._data.at[recipe.hash, self.key] = value
        else:
            return cached

        return value

    def histogram(
        self,
        progress: bool = False,
    ) -> "plotly.graph_objects.Figure":

        import plotly.graph_objects as go

        values = self.values

        fig = go.Figure(go.Histogram(x=values))

        fig.update_layout(xaxis_title=self.key, yaxis_title="count")

        return fig

    def unweighted(
        self,
        recipe: "Recipe",
    ) -> float:

        value = self.get_value(recipe)

        score = float(self.percentile_interpolator(value))

        if self.inverse:
            score = 1 - score

        return score

    ### DUNDERS

    def __call__(
        self,
        recipe: "Recipe",
    ) -> float:
        """return the score of a given value"""

        if not self.weight:
            return 0.0

        value = self.unweighted(recipe)

        # if value is None:
        # return 0.5

        return self.weight * value

    def __str__(self) -> str:
        """Unformatted string representation"""
        if self.weight is None:
            return f'{self._type}("{self.key}", inverse={self.inverse})'
        else:
            return f'{self._type}("{self.key}", weight={self.weight:.2f}, inverse={self.inverse})'

    def __repr__(self) -> str:
        """ANSI Formatted string representation"""
        import mcol

        return f"{mcol.bold}{mcol.underline}{self}{mcol.unbold}{mcol.ununderline}"

    def __rich__(self) -> str:
        """Rich Formatted string representation"""
        return f"[bold underline]{self}"


class CustomAttribute(Attribute):

    _type = "CustomAttribute"

    def __init__(self, scorer, key, function):
        self._function = function
        super(CustomAttribute, self).__init__(scorer=scorer, key=key)

    ### METHODS

    def get_value(
        self,
        recipe: "Recipe",
        serialise_price: bool = True,
        force: bool = False,
    ) -> float:

        if not force:
            cached = self.scorer._data[self.key][recipe.hash]

        if force or cached is None:

            value = self._function(recipe)

            if serialise_price and self.key == "price":
                value = value.amount
            self.scorer._data.at[recipe.hash, self.key] = value
        else:
            return cached

        return value


# DEFAULT_ATTRIBUTES = {
#     "num_scaffolds": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_compounds.count_by_tag(tag="Syndirella scaffold"),
#         description="The number of Syndirella scaffold compounds in this selection",
#     ),
#     "num_products": dict(
#         type="standard",
#         weight=1.0,
#         description="The number of product compounds in this selection",
#     ),
#     "num_scaffolds_elaborated": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_compounds.num_scaffolds_elaborated,
#         description="The number of Syndirella scaffold compounds that have at least one elaboration in this selection",
#     ),
#     "elaboration_balance": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_compounds.elaboration_balance,
#         description="A measure for how evenly scaffold compounds have been elaborated",
#     ),  ### REALLY UNPERFORMANT?
#     "num_inspirations": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_poses.num_inspirations,
#         description="The number of unique fragment compounds that inspired poses for product compounds in this selection",
#     ),
#     "num_inspiration_sets": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_poses.num_inspiration_sets,
#         description="The number of unique fragment combinations that inspired poses for product compounds in this selection",
#     ),
#     "risk_diversity": dict(
#         type="custom",
#         weight=0.0,
#         function=lambda r: r.product_compounds.risk_diversity,
#         description="A measure of how evenly spread the risk of elaborations are for each scaffold compound. Risk in this case refers to the number of atoms added",
#     ),
#     "interaction_count": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_interactions.num_features,
#         description="The number of protein features that are being interecated with in this selection",
#     ),
#     "interaction_balance": dict(
#         type="custom",
#         weight=0.0,
#         function=lambda r: r.product_interactions.per_feature_count_std,
#         description="A measure for how evenly protein features are being interacted with in this selection",
#     ),
#     "num_subsites": dict(
#         type="custom",
#         weight=1.0,
#         function=lambda r: r.product_poses.num_subsites,
#         description="Count the number of subsites that poses in this set come into contact with",
#     ),
#     "subsite_balance": dict(
#         type="custom",
#         weight=0.0,
#         function=lambda r: r.product_poses.subsite_balance,
#         description="Count the number of subsites that poses in this set come into contact with",
#     ),
#     "avg_distance_score": dict(
#         type="custom",
#         weight=-0.0,
#         function=lambda r: r.product_poses.avg_distance_score,
#         description="Average distance score (e.g. RMSD to fragment inspirations) for poses in this set",
#     ),
#     "avg_energy_score": dict(
#         type="custom",
#         weight=-0.0,
#         function=lambda r: r.product_poses.avg_energy_score,
#         description="Average energy score (e.g. binding ddG) for poses in this set",
#     ),
#     # "reaction_risk": dict(type='custom', weight=1.0, function=None),
#     # "pockets?": dict(type='custom', weight=1.0, function=None),
#     # "chemical_diversity": dict(type='custom', weight=1.0, function=None),
#     # "DMS/sequence_variability": dict(type='custom', weight=1.0, function=None),
# }

DEFAULT_ATTRIBUTES = {
    "num_scaffolds": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.combined_compounds.count_by_tag(tag="Syndirella scaffold"),
        description="The number of Syndirella scaffold compounds in this selection. Higher is better.",
    ),
    "num_compounds": dict(
        type="standard",
        weight=1.0,
        description="The number of product compounds in this selection. Higher is better.",
    ),
    "num_scaffolds_elaborated": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.combined_compounds.num_scaffolds_elaborated,
        description="The number of Syndirella scaffold compounds that have at least one elaboration in this selection. Higher is better.",
    ),
    "elaboration_balance": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.combined_compounds.elaboration_balance,
        description="A measure for how evenly scaffold compounds have been elaborated using an h-index. Higher is better.",
    ),  ### REALLY UNPERFORMANT?
    "num_inspirations": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.poses.num_inspirations,
        description="The number of unique fragment compounds that inspired poses for product compounds in this selection. Higher is better.",
    ),
    "num_inspiration_sets": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.poses.num_inspiration_sets,
        description="The number of unique fragment combinations that inspired poses for product compounds in this selection. Higher is better.",
    ),
    # "risk_diversity": dict(
    #     type="custom",
    #     weight=0.0,
    #     function=lambda r: r.combined_compounds.risk_diversity,
    #     description="A measure of how evenly spread the risk of elaborations are for each scaffold compound. Risk in this case refers to the number of atoms added. Higher is better",
    # ), # REMOVED BECAUSE IT DOES NOT NECESSARILY IMPROVE AS PRODUCTS ARE ADDED
    "interaction_count": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.interactions.num_features,
        description="The number of protein features that are being interecated with in this selection. Higher is better.",
    ),
    "interaction_balance": dict(
        type="custom",
        weight=0.0,
        function=lambda r: r.interactions.per_feature_count_hirsch,
        description="A measure for how evenly protein features are being interacted with in this selection using an h-index. Higher is better",
    ),
    "num_subsites": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.poses.num_subsites,
        description="Count the number of subsites that poses in this set come into contact with. Higher is better.",
    ),
    "subsite_balance": dict(
        type="custom",
        weight=0.0,
        function=lambda r: r.poses.subsite_balance,
        description="Count the number of subsites that poses in this set come into contact with",
    ),
    "avg_distance_score": dict(
        type="custom",
        weight=-0.0,
        function=lambda r: r.poses.avg_distance_score,
        description="Average distance score (e.g. RMSD to fragment inspirations) for poses in this set. Lower is better.",
    ),
    "avg_energy_score": dict(
        type="custom",
        weight=-0.0,
        function=lambda r: r.poses.avg_energy_score,
        description="Average energy score (e.g. binding ddG) for poses in this set. Lower is better.",
    ),
    # "reaction_risk": dict(type='custom', weight=1.0, function=None),
    # "pockets?": dict(type='custom', weight=1.0, function=None),
    # "chemical_diversity": dict(type='custom', weight=1.0, function=None),
    # "DMS/sequence_variability": dict(type='custom', weight=1.0, function=None),
}
