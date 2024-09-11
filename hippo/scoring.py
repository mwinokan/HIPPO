import mout
import numpy as np

# from .block import BuildingBlockSet
from scipy.interpolate import interp1d
import pandas as pd

import logging

logger = logging.getLogger("HIPPO")


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
        recipes: "RecipeSet",
        attributes: list[str],
    ) -> None:

        self._db = db
        self._recipes = recipes

        self._sorted_keys = None
        self._scores = None

        self._attributes = {}

        for key in attributes:
            attribute = Attribute(self, key)
            self._attributes[key] = attribute

        self.weights = 1.0

        self._df = None
        self._df_params = None

    ### FACTORIES

    @classmethod
    def default(
        cls,
        db: "Database",
        directory: "Path | str",
        pattern: str = "*.json",
    ):

        from .recipe import RecipeSet

        rset = RecipeSet(db, directory, pattern=pattern)

        self = cls.__new__(cls)

        attributes = [
            k for k, v in DEFAULT_ATTRIBUTES.items() if v["type"] == "standard"
        ]

        self.__init__(db=db, recipes=rset, attributes=attributes)

        # custom attributes
        for key, attribute in [
            (k, v) for k, v in DEFAULT_ATTRIBUTES.items() if v["type"] == "custom"
        ]:

            # print(key, attribute)

            self.add_custom_attribute(
                key, attribute["function"], weight_reset_warning=False
            )

        # weights

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

        self.__flag_modification()

        if isinstance(ws, float):
            ws = [ws] * self.num_attributes

        ws = [w for w in ws]
        wsum = sum([abs(w) for w in ws])

        for a, w in zip(self.attributes, ws):
            a.weight = w / wsum

    @property
    def scores(self):
        if self._scores is None:
            logger.debug("Scorer.scores")
            self._scores = {k: self.score(r) for k, r in self.recipes.items()}

        return self._scores

    @property
    def sorted_keys(self):

        if self._sorted_keys is None:
            logger.debug("Scorer.sorted_keys")
            self._sorted_keys = [
                t[0]
                for t in sorted(
                    [t for t in self.scores.items()], key=lambda t: t[0], reverse=True
                )
            ]

        return self._sorted_keys

    @property
    def sorted(self):
        return [self.recipes[key] for key in self.sorted_keys]

    @property
    def best(self):
        key = self.sorted_keys[0]
        return self.recipes[key]

    ### METHODS

    def add_custom_attribute(
        self,
        key,
        function,
        weight_reset_warning: bool = True,
    ) -> "CustomAttribute":

        ca = CustomAttribute(self, key, function)

        if key not in self._attributes:

            self.__flag_modification()

            self._attributes[key] = ca

            if weight_reset_warning:
                logger.warning("Attribute weights have been reset")
            self.weights = 1.0

        else:
            logger.warning("Existing attribute with {key=}")

        return self._attributes[key]

    def score(
        self,
        recipe: "Recipe",
        *,
        debug: bool = False,
    ) -> float:

        score = sum([a(recipe) for a in self.attributes])

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
            logger.var("score", score)

        recipe._score = score

        return score

    def get_df(
        self,
        serialise_price: bool = True,
        debug: bool = True,
        **kwargs,
    ) -> "pandas.DataFrame":

        from pandas import DataFrame

        params = dict(serialise_price=serialise_price, **kwargs)

        if self._df is None or params != self._df_params:

            if debug:
                logger.debug("Scorer.get_df()")

            if debug:
                logger.debug("Scorer.recipes.get_df()")

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
                    logger.debug(f"Getting {attribute.key=} values")

                if isinstance(attribute, CustomAttribute):

                    df[attribute.key] = attribute.values

                else:

                    df[attribute.key] = self.recipes.get_values(
                        key=attribute.key, serialise_price=serialise_price
                    )

            self._df = df
            self._df_params = params

        return self._df

    def plot(
        self,
        keys: str | list[str],
        **kwargs,
    ):

        import plotly.express as px

        df = self.get_df(**kwargs)

        if isinstance(keys, str):
            return px.histogram(df, x=keys)

        if len(keys) != 2:
            logger.error("Only two keys supported")
            return None

        return px.scatter(df, x=keys[0], y=keys[1], color="score")

    def top(self, n):
        return [self.recipes[key] for key in self.sorted_keys[:n]]

    ### INTERNALS

    def __flag_modification(self):
        self._df = None
        self._df_params = None
        self._scores = None
        self._sorted_keys = None

    def summary(self):

        logger.header(repr(self))
        for attribute in self.attributes:
            logger.out(
                f"{repr(attribute)} min={attribute.min:.3g}, mean={attribute.mean:.3g}, std={attribute.std:.3g}, max={attribute.max:.3g}"
            )

    ### DUNDERS

    def __repr__(self):
        return f"Scorer(#recipes={self.num_recipes})"


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

        self._mean = None
        self._std = None
        self._min = None
        self._max = None

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
        if not self._value_dict:
            logger.debug(f"Attribute(key={self.key}).value_dict")
            for recipe in self.scorer.recipes:
                self._value_dict[recipe.hash] = self.get_value(recipe)
        return self._value_dict

    @property
    def values(self):
        return list(self.value_dict.values())

    @property
    def mean(self):
        if self._mean is None:
            self._mean = np.mean(self.values)
        return self._mean

    @property
    def std(self):
        if self._std is None:
            self._std = np.std(self.values)
        return self._std

    @property
    def max(self):
        if self._max is None:
            self._max = max(self.values)
        return self._max

    @property
    def min(self):
        if self._min is None:
            self._min = min(self.values)
        return self._min

    @property
    def weight(self):
        return self._weight

    @weight.setter
    def weight(self, w):
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

    def __repr__(self) -> str:

        if self.weight is None:
            return f'{self._type}("{self.key}", inverse={self.inverse})'
        else:
            return f'{self._type}("{self.key}", weight={self.weight:.2f}, inverse={self.inverse})'

    ### METHODS

    def get_value(
        self,
        recipe: "Recipe",
        serialise_price: bool = True,
    ) -> float:
        value = getattr(recipe, self.key)
        if serialise_price and self.key == "price":
            value = value.amount
        return value

    def histogram(
        self,
        progress: bool = False,
    ) -> "plotly.graph_objects.Figure":

        import plotly.graph_objects as go

        values = self.scorer.recipes.get_values(
            self.key, progress=progress, serialise_price=True
        )

        fig = go.Figure(go.Histogram(x=values))

        fig.update_layout(xaxis_title=self.key, yaxis_title="count")

        return fig

    def unweighted(
        self,
        recipe: "Recipe",
        value: float = None,
        debug: bool = False,
    ) -> float:

        if value is None:
            value = self.get_value(recipe)

        if debug:
            logger.debug(f"{value=}")

        score = float(self.percentile_interpolator(value))

        if debug:
            logger.debug(f"{score=}")

        if self.inverse:
            score = 1 - score

        return score


class CustomAttribute(Attribute):

    _type = "CustomAttribute"

    def __init__(self, scorer, key, function):
        self.get_value = function
        self._values = None
        super(CustomAttribute, self).__init__(scorer=scorer, key=key)

    ### PROPERTIES

    @property
    def values(self):

        if self._values is None:

            from tqdm import tqdm

            logger.debug(f"CustomAttribute(key={self.key}).values")

            self._values = []
            for recipe in tqdm(self.scorer.recipes):
                value = self.get_value(recipe)
                self._values.append(value)

        return self._values

    ### METHODS

    def histogram(
        self,
        progress: bool = False,
    ) -> "plotly.graph_objects.Figure":

        import plotly.graph_objects as go

        fig = go.Figure(go.Histogram(x=self.values))

        fig.update_layout(xaxis_title=self.key, yaxis_title="count")

        return fig


DEFAULT_ATTRIBUTES = {
    "num_bases": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_compounds.count_by_tag(tag="Syndirella base"),
        description="The number of Syndirella base compounds in this selection",
    ),
    "num_products": dict(
        type="standard",
        weight=1.0,
        description="The number of product compounds in this selection",
    ),
    "num_bases_elaborated": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_compounds.num_bases_elaborated,
        description="The number of Syndirella base compounds that have at least one elaboration in this selection",
    ),
    "elaboration_balance": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_compounds.elaboration_balance,
        description="A measure for how evenly base compounds have been elaborated",
    ),
    "num_inspirations": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_poses.num_inspirations,
        description="The number of unique fragment compounds that inspired poses for product compounds in this selection",
    ),
    "num_inspiration_sets": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_poses.num_inspiration_sets,
        description="The number of unique fragment combinations that inspired poses for product compounds in this selection",
    ),
    "risk_diversity": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_compounds.risk_diversity,
        description="A measure of how evenly spread the risk of elaborations are for each base compound. Risk in this case refers to the number of atoms added",
    ),
    "interaction_count": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_interactions.num_features,
        description="The number of protein features that are being interecated with in this selection",
    ),
    "interaction_balance": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_interactions.per_feature_count_std,
        description="A measure for how evenly protein features are being interacted with in this selection",
    ),
    "num_subsites": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_poses.num_subsites,
        description="Count the number of subsites that poses in this set come into contact with",
    ),
    "subsite_balance": dict(
        type="custom",
        weight=1.0,
        function=lambda r: r.product_poses.subsite_balance,
        description="Count the number of subsites that poses in this set come into contact with",
    ),
    # "reaction_risk": dict(type='custom', weight=1.0, function=None),
    # "pockets?": dict(type='custom', weight=1.0, function=None),
    # "chemical_diversity": dict(type='custom', weight=1.0, function=None),
    # "DMS/sequence_variability": dict(type='custom', weight=1.0, function=None),
}
