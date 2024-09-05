
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

    def __init__(self, 
        db: 'Database', 
        recipes: 'RecipeSet', 
        attributes: list[str], 
    ) -> None:

        self._db = db
        self._recipes = recipes

        # self._sorted = None

        self._attributes = []

        for key in attributes:
            attribute = Attribute(self, key)
            self._attributes.append(attribute)

        self.weights = 1.0

    ### PROPERTIES

    @property
    def num_recipes(self) -> int:
        """Number of recipes being evaluated"""
        return len(self._recipes)

    @property
    def attributes(self):
        return self._attributes

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

        if isinstance(ws, float):
            ws = [ws] * self.num_attributes

        ws = [w for w in ws]
        wsum = sum([abs(w) for w in ws])

        for a, w in zip(self.attributes, ws):
            a.weight = w / wsum

    @property
    def scores(self):
        return [self.score(r) for r in self.recipes]

#     @property
#     def sorted(self):
#         if self.verbosity:
#             mout.debug(f"sorting {len(self.bb_sets)} BB sets...")

#         if self._sorted is None:
#             self._sorted = sorted(
#                 self.bb_sets, key=lambda x: self.score(x), reverse=True
#             )
#         return self._sorted

#     @property
#     def best(self):
#         return self.sorted[0]

#     ### METHODS

    # @property
    # def get_values(self):
    #     return self._get_values


#     def add_attribute(self, a):
#         self._attributes.append(a)

    def score(self, 
        recipe: 'Recipe', 
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

        return score

    def get_df(self,
        serialise_price: bool = True,
        **kwargs,
    ) -> 'pandas.DataFrame': 

        logger.debug('Scorer.get_df()')

        df = self.recipes.get_df(serialise_price=serialise_price, **kwargs)

        df['score'] = self.scores

        for attribute in self.attributes:
            df[attribute.key] = self.recipes.get_values(key=attribute.key, serialise_price=serialise_price)

        return df

    def plot(self, 
        keys: list[str],
        **kwargs,
    ):

        import plotly.express as px

        return px.scatter(self.get_df(**kwargs), x=keys[0], y=keys[1], color='score')

#     def top(self, n):
#         return self.sorted[:n]

    ### DUNDERS

    def __repr__(self):
        return f'Scorer(#recipes={self.num_recipes})'

class Attribute:

    _type = 'Attribute'

    ### DUNDERS

    def __init__(self, 
        scorer: 'Scorer',
        key: str, 
        *,
        inverse: bool = False, 
        weight: float = 1.0,
        bins: int = 100,
    ):

        self._scorer = scorer
        
        self._key = key
        self._inverse = None
        self._weight = weight
        
        self._value_dict = {}

        self._mean = None
        self._std = None

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
            for recipe in self.scorer.recipes:
                self._value_dict[recipe.hash] = self.get_value(recipe)
        return self._value_dict

    @property
    def values(self):
        return list(self.value_dict.values())

    @property
    def mean(self):
        if self._mean is None:
            self._mean = np.mean(self.values.values())
        return self._mean

    @property
    def std(self):
        if self._std is None:
            self._std = np.std(self.values.values())
        return self._std

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

    def __call__(self, 
        recipe: 'Recipe',
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

    def get_value(self, 
        recipe: 'Recipe',
        serialise_price: bool = True,
    ) -> float:
        value = getattr(recipe, self.key)
        if serialise_price and self.key == 'price':
            value = value.amount
        return value

    def histogram(self, 
        progress: bool = False,
    ) -> 'plotly.graph_objects.Figure':

        import plotly.graph_objects as go

        values = self.scorer.recipes.get_values(self.key, progress=progress, serialise_price=True)

        fig = go.Figure(go.Histogram(x=values))

        fig.update_layout(xaxis_title=self.key, yaxis_title='count')

        return fig

    def unweighted(self, 
        recipe: 'Recipe', 
        value: float = None,
        debug: bool = False,
    ) -> float:

        if value is None:
            value = self.get_value(recipe)

        if debug:
            logger.debug(f'{value=}')

        score = float(self.percentile_interpolator(value))
        
        if debug:
            logger.debug(f'{score=}')

        if self.inverse:
            score = 1 - score

        return score


class CustomAttribute(Attribute):

    # _type = 'CustomAttribute'

    def __init__(self, key, bb_sets, function):
        raise NotImplementedError
#         self.get_value = function
#         super(CustomAttribute, self).__init__(key=key, bb_sets=bb_sets)
#         self._type = "CustomAttribute"
