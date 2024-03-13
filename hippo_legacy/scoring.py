
import mout
import numpy as np
from .block import BuildingBlockSet
from scipy.interpolate import interp1d
import pandas as pd

class Attribute:

    ### DUNDERS
    
    def __init__(self, key, bb_sets, bins=100):
        self.key = key
        self.reverse = None
        self.weight = None
        self._type = 'Attribute'

        values = [self.get_value(bbs) for bbs in bb_sets]

        self.mean = np.mean(values)
        self.std = np.std(values)
        
        self.get_percentile_interpolator(values,bins=bins)
        
    def __call__(self,value):
        """return the score of a given value"""

        if not self.weight:
            return 0.0

        if isinstance(value,BuildingBlockSet) or 'BuildingBlockSet' in str(type(value)):
            value = self.get_value(value)
        elif isinstance(value,dict):
            value = self.get_value(value)

        value = float(value)

        if value is None:
            return 0.5

        return self.weight * self.unweighted(value)

    def __repr__(self):
        return f'{self._type}("{self.key}", weight={self.weight:.2f}, reverse={self.reverse})'

    ### METHODS

    def get_value(self,bb_set):
        if isinstance(bb_set,dict):
            return bb_set[self.key]
        else:
            return getattr(bb_set,self.key)

    def get_percentile_interpolator(self,values,bins=100):
        count, bins_count = np.histogram(values, bins=bins)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        self._percentile_interpolator = interp1d(bins_count[1:],cdf,kind='linear',fill_value="extrapolate")

    def unweighted(self,value):
        try:
            if self.reverse:
                return 1 - self._percentile_interpolator(value)
            else:
                return self._percentile_interpolator(value)
        except ValueError as e:
            mout.error(f'{self.key=}')
            mout.error(f'{value=}')
            mout.error(f'{type(value)=}')
            raise e

class CustomAttribute(Attribute):

    def __init__(self, key, bb_sets, function):
        self.get_value = function
        super(CustomAttribute,self).__init__(key=key, bb_sets=bb_sets)
        self._type = 'CustomAttribute'

class Scorer:

    """

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
    
    def __init__(self, animal, bb_sets, attributes, verbosity=1):

        if verbosity:
            mout.debug(f'Scorer({animal})')

        self.animal = animal
        self.bb_sets = bb_sets
        self.verbosity = verbosity

        self._sorted = None

        self._attributes = []

        for key in attributes:
            attribute = Attribute(key, bb_sets)
            self._attributes.append(attribute)

    ### PROPERTIES

    @property
    def attributes(self):
        return self._attributes

    @property
    def num_attributes(self):
        return len(self.attributes)

    @property
    def weights(self):
        return [a.weight for a in self.attributes]

    @weights.setter
    def weights(self, ws):
        
        rs = [w < 0 for w in ws]
        ws = [abs(w) for w in ws]

        for a,w,r in zip(self.attributes,ws,rs):
            a.weight = w / sum(ws)
            a.reverse = r

    @property
    def sorted(self):
        if self.verbosity:
            mout.debug(f'sorting {len(self.bb_sets)} BB sets...')

        if self._sorted is None:
            self._sorted = sorted(self.bb_sets, key=lambda x: self.score(x), reverse=True)
        return self._sorted

    @property
    def best(self):
        return self.sorted[0]
    
    ### METHODS

    def add_attribute(self,a):
        self._attributes.append(a)

    def score(self, bb_set, verbosity=0):

        score = sum([a(bb_set) for a in self.attributes])

        if verbosity:
            print_data = []
            for attribute in self.attributes:
                print_data.append(dict(
                    key=attribute.key,
                    weight=attribute.weight,
                    value=attribute.get_value(bb_set),
                    unweighted=attribute.unweighted(attribute.get_value(bb_set)),
                    weighted=attribute(bb_set)
                ))

            print(pd.DataFrame(print_data))
            mout.var('score',score)

        if isinstance(bb_set,dict):
            bb_set['score'] = score
        else:
            bb_set.score = score

        return score

    def top(self, n):
        return self.sorted[:n]
