
from .set import CompoundSet
import pandas as pd
import mout
import json
import pickle
import numpy as np
from collections.abc import MutableSet
from pprint import pprint
from rdkit.Chem import CanonSmiles
import plotly.graph_objects as go
from scipy.interpolate import interp1d
from copy import deepcopy
import molparse as mp
from pprint import pprint
import functools

class PriceInterpolator:

    def __init__(self,log_price_interpolator, min_amount, min_price):
        self.log_price_interpolator = log_price_interpolator
        self.min_amount = min_amount
        self.min_price = min_price

    def __call__(self,x):

        if isinstance(x,np.ndarray):
            return [ self.__call__(v) for v in x]

        if x <= self.min_amount:
                return self.min_price
        else:
            return np.power(10,self.log_price_interpolator(np.log10(x)))

class FlatPriceInterpolator:

    def __init__(self, unit_price, unit_in_mg):
        self.price_per_mg = unit_price / unit_in_mg
        self.min_amount = unit_in_mg
        self.min_price = unit_price

    def __call__(self,x):
        if isinstance(x,np.ndarray):
            return [ self.__call__(v) for v in x]

        if x <= self.min_amount:
            return self.min_price
        else:
            return self.price_per_mg * x

class SingleValueInterpolator:

    def __init__(self, value):
        self.value = value

    def __call__(self,x):
        if isinstance(x,np.ndarray):
            return [ self.__call__(v) for v in x]

        return self.value

class BuildingBlock:

    ### DUNDERS

    def __init__(self, smiles, molecular_weight=None):
        
        self._smiles = smiles
        self._molecular_weight = molecular_weight
        self._purchase_info = []
        self._has_purchase_interpolators = False
        self._min_lead_time = None
        self._required_amount = None
        self._price_interpolator = None
        self._has_purchase_interpolators = False
        self._lead_time_interpolator = None

    def __repr__(self):
        # if self._cost_str is not None:
        #     return f'BuildingBlock({self.smiles}, cost={self._cost_str}, lead_time={self.lead_time} weeks)'
        # else:
        return f'BuildingBlock({self.smiles})'

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __hash__(self):
        return hash(self.smiles)

    def __deepcopy__(self):
        copy = BuildingBlock(self.smiles, self._molecular_weight)
        copy._purchase_info = deepcopy(self._purchase_info)
        copy._has_purchase_interpolators = self._has_purchase_interpolators
        copy._min_lead_time = self._min_lead_time
        copy._required_amount = self._required_amount
        copy._price_interpolator = self._price_interpolator
        copy._has_purchase_interpolators = self._has_purchase_interpolators
        copy._lead_time_interpolator = self._lead_time_interpolator
        # copy._products = self._products
        return copy

    ### PROPERTIES

    @property
    def smiles(self):
        return self._smiles

    @property
    def purchaseable(self):
        return bool(self._purchase_info)

    @property
    def has_purchase_interpolators(self):
        return self._has_purchase_interpolators

    @property
    def min_lead_time(self):
        return self._min_lead_time

    @property
    def required_amount(self):
        return self._required_amount
    
    @required_amount.setter
    def required_amount(self, a):
        self._required_amount = float(a)
    
    ### METHODS

    def copy(self):
        return self.__deepcopy__()

    @functools.cache
    def get_purchase_info(self,amount_in_mg):

        for info in self._purchase_info:
            if info['amount'] == amount_in_mg:
                return info
        
        if self.has_purchase_interpolators:
            return dict(price=self._price_interpolator(amount_in_mg),lead_time=self._lead_time_interpolator(amount_in_mg))
        else:
            mout.error(f'No purchase info for {amount_in_mg}mg of {self}')
            return {}

    @functools.cache
    def get_price(self,amount_in_mg=10):
        if amount_in_mg == 0:
            return 0
        info = self.get_purchase_info(amount_in_mg)
        if not info:
            return None
        return info['price']

    @functools.cache
    def get_lead_time(self,amount_in_mg=10):
        info = self.get_purchase_info(amount_in_mg)
        return float(info['lead_time'])

    def add_purchase_info(self,**info):

        assert 'lead_time' in info
        assert 'price' in info
        assert 'amount' in info
        assert info['amount_unit'] == 'mg'
        assert info['price_unit'] == 'USD'
        assert info['lead_time_unit'] == 'days'

        self._purchase_info.append(info)

    def create_price_interpolator(self, log_price_interpolator, min_amount, min_price):

        def price_interpolator(x):
            if x < min_amount:
                return min_price
            else:
                return np.power(10,log_price_interpolator(np.log10(x)))

        return price_interpolator

    def generate_purchase_interpolators(self, debug=False, show=False):
        
        self._purchase_info = sorted(self._purchase_info, key=lambda x: x['amount'])

        amounts = []
        prices = []
        lead_times = []

        for info in self._purchase_info:

            price = info['price']
            amount = info['amount']
            lead_time = info['lead_time']
            
            amounts.append(amount)
            lead_times.append(lead_time)
            prices.append(price)

        if debug:
            mout.header(self)
            mout.var('amounts',amounts)
            mout.var('lead_times',lead_times)
            mout.var('prices',prices)

        if len(amounts) < 2:
            mout.warning(f'Single pricing data point ({amount}mg) for {self.smiles}')
            self._lead_time_interpolator = SingleValueInterpolator(lead_time)
            self._price_interpolator = FlatPriceInterpolator(price,amount)

        else:

            log_amounts = np.log10(amounts)
            log_prices = np.log10(prices)

            log_price_interpolator = interp1d(log_amounts,log_prices,kind='linear',fill_value='extrapolate')

            self._price_interpolator = PriceInterpolator(log_price_interpolator, amounts[0], prices[0])

            self._lead_time_interpolator = interp1d(amounts,lead_times,kind='linear',bounds_error=False,fill_value=(lead_times[0],lead_times[-1]))

            self._min_lead_time = min(lead_times)

        if show:

            fig = go.Figure()
            
            trace = go.Scatter(name='price',x=amounts,y=prices,mode='markers')
            fig.add_trace(trace)

            trace = go.Scatter(name='lead_time',x=amounts,y=lead_times,mode='markers')
            fig.add_trace(trace)

            fig.update_xaxes(type="log")
            fig.update_yaxes(type="log")
                
            # interpolators
            xs = np.linspace(0.1,2000,500)
            
            ys = self._price_interpolator(xs)
            trace = go.Scatter(name='price (interpolated)',x=xs,y=ys,mode='lines')
            fig.add_trace(trace)

            ys = self._lead_time_interpolator(xs)
            trace = go.Scatter(name='lead_time (interpolated)',x=xs,y=ys,mode='lines')
            fig.add_trace(trace)

            fig.show()

        self._has_purchase_interpolators = True

class BuildingBlockSet(MutableSet):

    ### DUNDERS

    def __init__(self, bbs=()):
        # super(BuildingBlockSet,self).__init__(s)
        self._name = None
        self._cost = None

        self._elements = []

        for bb in bbs:
            if bb not in self._elements:
                self._elements.append(bb)

        self._products = None

    def __iter__(self):
        return iter(self._elements)

    def __contains__(self, query):
        if isinstance(query,str):
            return query in self.smiles
        elif isinstance(query,BuildingBlock):
            return query in self._elements
        else:
            return NotImplementedError(f"Unsupported query type: {type(query)}")

    def __len__(self):
        return len(self._elements)

    def __getitem__(self, key):

        if isinstance(key,list):
            return BuildingBlockSet([bb for bb in [self[k] for k in key] if bb is not None])

        if isinstance(key,BuildingBlock):
            key = key.smiles

        key = CanonSmiles(key)

        matches = [bb for bb in self if bb.smiles == key]

        if len(matches) < 1:
            mout.error(f'{key} not in {self}')
            return None
        elif len(matches) > 1:
            mout.error(f'Multiple {key} in {self}')
            exit()
            return None

        return matches[0]

    def __repr__(self):
        if self.name is not None:
            return f'BuildingBlockSet("{self.name}", #building_blocks={len(self)})'
        else:
            return f'BuildingBlockSet(#{self.id}, #building_blocks={len(self)})'

    def __hash__(self):
        return hash(" ".join(sorted(self.smiles)))

    ### PROPERTIES

    @property
    def id(self):
        return hash(self)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    # @property
    # def cost(self):
    #     return self._cost

    @property
    def smiles(self):
        return[ bb.smiles for bb in self._elements]

    @property
    def amounts(self):
        return {bb.smiles:bb.required_amount for bb in bb_set}

    @property
    def products(self):
        if self._products is None:
            mout.error('Products have not been generated')
        return self._products

    @property
    def num_products(self):
        return len(self.products)

    @property
    def hit_feature_coverage(self):
        return self._hit_feature_coverage

    @property
    def num_new_features(self):
        return self._num_new_features
    
    ### METHODS

    def copy(self):
        new = BuildingBlockSet([bb.copy() for bb in self])
        return new

    def discard(self, key):
        if key in self:
            i = self._elements.index(key)
            del self._elements[i]
        else:
            mout.warning(f'Tried to delete element not in set: {key}')

    def add(self, bb):
        if bb not in self._elements:
            self._elements.append(bb)
            bb._set_name = self.name

    def get_products(self, candidates):

        comp_set = CompoundSet(f'Products({self})')

        candidates = [c for c in candidates if c.building_blocks is not None]

        for candidate in candidates:

            if any([bb not in self for bb in candidate.building_blocks]):
                continue

            comp_set.add(candidate)

        comp_set.immutable = True

        # get required amounts
        for bb in self:
            bb.required_amount = 0
        for prod in comp_set:
            for bb in prod.building_blocks:
                if bb in self:
                    self[bb].required_amount += 1

        self._products = comp_set

        return comp_set

    def summary(self):

        mout.header(self)

        print_data = []
        for bb in self:
            print_data.append(dict(
                smiles=bb.smiles,
                cost=bb.get_price(),
                lead_time=bb.get_lead_time(),
            ))

        print(pd.DataFrame(print_data))

    def assign_price_info(self,minimise='cost_min'):

        for bb in self:

            if bb._catalog_metadata is not None:            

                if minimise == 'cost_min':
                    bb.get_lowest_cost()
                elif minimise == 'cost_max':
                    bb.get_lowest_cost(use_max=True)
                elif minimise == 'lead_time':
                    bb.get_shortest_lead_time()
                else:
                    raise Exception(f"Unknown quantity to minimise '{minimise}'")

    def get_price(self,blanket_amount=None):

        price = 0.0

        for bb in self._elements:
            
            if blanket_amount is not None:
                amount = blanket_amount
            else:
                amount = bb.required_amount

            assert not np.isnan(amount), (bb, amount)

            bb_price = bb.get_price(amount)

            if bb_price is not None:

                assert not np.isnan(bb_price), (bb, amount, bb_price, bb.generate_purchase_interpolators(debug=True,show=True))

                price += bb_price

        return price

    def get_lead_times(self,blanket_amount=None):

        lead_times = []

        for bb in self._elements:
            
            if blanket_amount:
                amount = blanket_amount
            else:
                amount = bb.required_amount

            lead_times.append(bb.get_lead_time(amount))

        return lead_times

    def remove_unused(self,candidates,verbosity=0):
        self.get_products(candidates)

        count = 0
        while len([bb for bb in self if bb.required_amount < 1]) > 0:
            for bb in self:
                if bb.required_amount < 1:
                    self.remove(bb)
                    if verbosity:
                        mout.warning(f'Removing unused {bb}')
                    count += 1
                self.get_products(candidates)
        return count

    def write(self,directory):

        identifier = self.id

        mp.write(directory / f'BBS_{identifier}.pickle',self)

        json_data = dict(id=identifier,size=len(self),price=self.get_price(),building_blocks=[],products=[])

        for bb in self:
            json_data['building_blocks'].append(dict(smiles=bb.smiles, amount=bb.required_amount, price=bb.get_price(bb.required_amount), lead_time=bb.get_lead_time(bb.required_amount)))

        try:
            getattr(self, '_products')
        except AttributeError:
            self.get_products()

        for prod in self._products:
            json_data['products'].append(dict(name=prod.name, smiles=str(prod.smiles)))

        json_data['num_products'] = len(json_data['products'])

        mp.write(directory / f'BBS_{identifier}.json', json_data)

    def prepare_for_scoring(self, hit_features):

        self_features = self.products.get_present_features()

        num_in_hits_but_not_self = len(hit_features - self_features)
        self._hit_feature_coverage = (len(hit_features)-num_in_hits_but_not_self)/len(hit_features)

        self._num_new_features = len(self_features - hit_features)
