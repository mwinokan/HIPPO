
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

from .compound import FailedToAssignBondOrders

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
        self._molecular_weight = float(molecular_weight)
        self._purchase_info = []
        self._has_purchase_interpolators = False
        self._min_lead_time = None
        self._price_1mg = None
        self._price_10mg = None
        self._lead_time_1mg = None
        self._lead_time_10mg = None
        self._required_amount = None
        self._price_interpolator = None
        self._has_purchase_interpolators = False
        self._lead_time_interpolator = None
        self._id = None

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
        copy._id = self._id
        copy._purchase_info = deepcopy(self._purchase_info)
        copy._min_lead_time = self._min_lead_time
        copy._required_amount = self._required_amount
        copy._price_interpolator = self._price_interpolator
        copy._has_purchase_interpolators = self._has_purchase_interpolators
        copy._lead_time_interpolator = self._lead_time_interpolator
        copy._price_1mg = self._price_1mg
        copy._price_10mg = self._price_10mg
        copy._lead_time_1mg = self._lead_time_1mg
        copy._lead_time_10mg = self._lead_time_10mg

        # copy._products = self._products
        return copy

    ### PROPERTIES

    @property
    def smiles(self):
        return self._smiles

    @property
    def molecular_weight(self):
        return self._molecular_weight

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self,i):
        self._id = i

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
    def price_1mg(self):
        return self._price_1mg

    @property
    def price_10mg(self):
        return self._price_10mg

    @property
    def lead_time_1mg(self):
        return self._lead_time_1mg

    @property
    def lead_time_10mg(self):
        return self._lead_time_10mg

    @property
    def required_amount(self):
        return self._required_amount
    
    @property
    def dict(self):
        return dict(
            id = self.id,
            smiles = self.smiles,
            molecular_weight = self.molecular_weight,
            purchaseable = self.purchaseable,
            price_1mg = self.price_1mg,
            price_10mg = self.price_10mg,
            lead_time_1mg = self.lead_time_1mg,
            lead_time_10mg = self.lead_time_10mg,
            min_lead_time = self.min_lead_time,
            required_amount = self.required_amount,
        )
    
    ### METHODS

    def clear_amount(self):
        self._required_amount = 0

    def increment_amount(self):
        self._required_amount += 1

    def copy(self):
        return self.__deepcopy__()

    # @functools.cache
    def get_purchase_info(self,amount_in_mg,verbosity=1):

        for info in self._purchase_info:
            if info['amount'] == amount_in_mg:
                return info
        
        if self.has_purchase_interpolators:
            return dict(price=self._price_interpolator(amount_in_mg),lead_time=self._lead_time_interpolator(amount_in_mg))
        else:
            if verbosity:
                mout.error(f'No purchase info for {amount_in_mg}mg of {self}')
            return {}

    def replace_purchase_info(self,amount_in_mg,new):

        for i,info in enumerate(self._purchase_info):
            
            if info['amount'] == amount_in_mg:
                break

        else:
            mout.error(f'could not replace purchase info {amount_in_mg}mg: {self}')
            return

        self._purchase_info[i] = new
        # self._purchase_info.append(new)

    # @functools.cache
    def get_price(self,amount_in_mg=10):
        if amount_in_mg == 0:
            return 0
        info = self.get_purchase_info(amount_in_mg)
        if not info:
            return None
        return info['price']

    # @functools.cache
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

        amount = info['amount']

        amounts = [ i['amount'] for i in self._purchase_info]

        if amount in amounts:

            existing = self.get_purchase_info(amount)

            if all([info[key] == existing[key] for key in ['lead_time','price']]):
                return

            if all([info[key] >= existing[key] for key in ['lead_time','price']]):
                return

            if all([info[key] <= existing[key] for key in ['lead_time','price']]):
                self.replace_purchase_info(amount,info)
                return

            if info['price'] < existing['price']:
                self.replace_purchase_info(amount,info)
                mout.warning(f'Preferring ${info["price"]} & {info["lead_time"]}days to ${existing["price"]} & {existing["lead_time"]}days ({self.smiles})')
                return                
            else:
                mout.warning(f'Preferring ${existing["price"]} & {existing["lead_time"]}days to ${info["price"]} & {info["lead_time"]}days ({self.smiles})')
                return

            mout.error(f'Could not solve duplicate purchase info for {info["amount"]}mg: {self}')
            return

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

        self._price_1mg = self.get_price(1)
        self._price_10mg = self.get_price(10)
        self._lead_time_1mg = self.get_lead_time(1)
        self._lead_time_10mg = self.get_lead_time(10)

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
        self._hit_feature_coverage = None
        self._num_new_features = None
        
        self._avg_fragmenstein_ddG = None
        self._std_fragmenstein_ddG = None
        self._avg_fragmenstein_mRMSD = None
        self._std_fragmenstein_mRMSD = None

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
            # mout.error(f'{key} not in {self}')
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
            return f'BuildingBlockSet(#building_blocks={len(self)})'

    def __hash__(self):
        return hash(" ".join(sorted(self.smiles)))

    ### PROPERTIES

    @property
    def hash(self):
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
        return {bb.smiles:bb.required_amount for bb in self}

    @property
    def total_bb_amount(self):
        return sum(bb.required_amount for bb in self)

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
    
    @property
    def avg_fragmenstein_ddG(self):
        return self._avg_fragmenstein_ddG

    @property
    def std_fragmenstein_ddG(self):
        return self._std_fragmenstein_ddG

    @property
    def avg_fragmenstein_mRMSD(self):
        return self._avg_fragmenstein_mRMSD

    @property
    def std_fragmenstein_mRMSD(self):
        return self._std_fragmenstein_mRMSD

    @property
    def price(self):
        return self.get_price()

    @property
    def lead_times(self):
        return self.get_lead_times()

    @property
    def min_lead_time(self):
        return min(self.lead_times)

    @property
    def max_lead_time(self):
        return max(self.lead_times)

    @property
    def avg_lead_time(self):
        return np.mean(self.lead_times)

    @property
    def std_lead_time(self):
        return np.std(self.lead_times)

    @property
    def dict(self):

        d = dict(

            # properties
            # id = self.id,
            name = self.name,
            
            num_bbs=len(self),
            total_bb_amount=self.total_bb_amount,
            num_products=self.num_products,
            
            price=self.price,
            max_lead_time=self.max_lead_time,
            avg_lead_time=self.avg_lead_time,

            hit_feature_coverage=self.hit_feature_coverage,
            num_new_features=self.num_new_features,
            
            avg_fragmenstein_ddG=self.avg_fragmenstein_ddG,
            std_fragmenstein_ddG=self.std_fragmenstein_ddG,
            avg_fragmenstein_mRMSD=self.avg_fragmenstein_mRMSD,
            std_fragmenstein_mRMSD=self.std_fragmenstein_mRMSD,

            # lists
            building_blocks = [bb.dict for bb in self],
            products = [p.dict for p in self.products],
            
        )

        return d
    
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

    def get_products(self, candidates, debug=False):

        comp_set = CompoundSet(f'Products({self})')

        candidates = [c for c in candidates if c.building_blocks is not None]

        for candidate in candidates:

            if any([bb not in self for bb in candidate.building_blocks]):
                continue

            comp_set.add(candidate)

        comp_set.immutable = True

        # get required amounts
        for bb in self:
            bb.clear_amount()

        if debug:
            mout.var('#products',len(comp_set))

        count = 0
        for prod in comp_set:

            if debug:
                mout.out(f'BBS has all BBs for {prod.name}')

            for bb in prod.building_blocks:
                assert bb in self
                self[bb].increment_amount()
                count += 1

        self._products = comp_set

        return comp_set

    def summary(self,return_df=False):

        mout.header(self)

        print_data = []
        for bb in self:
            d = dict(
                smiles=bb.smiles,
                amount=bb.required_amount,
                # cost_per_mg=round(bb.get_price(),2),
            )
            
            d[f'price'] = round(bb.get_price(bb.required_amount),2)
            d[f'lead_time'] = round(bb.get_lead_time(bb.required_amount),1)

            print_data.append(d)

        df = pd.DataFrame(print_data)

        if return_df:
            return df

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

    def write(self,directory,bbs_id=None,extra_data=None,verbosity=1):

        identifier = bbs_id if bbs_id is not None else self.id

        self.name = f'BBS_{identifier:06}'

        mp.write(directory / f'BBS_{identifier:06}.pickle',self,verbosity=verbosity-1)

        json_data = self.dict

        if extra_data:
            json_data.update(extra_data)

        mp.write(directory / f'BBS_{identifier:06}.json', json_data,verbosity=verbosity-1)

    def product_df(self, animal):

        def get_hit(name):
            name = name.replace('-','_')
            for hit in animal.hits:
                if name in hit._pose_name:
                    return hit

        data = []

        for p in self.products:

            d = p.dict

            d['mol'] = p.mol

            for i,h in enumerate(p.inspirations):
                hit = get_hit(h)
                try:
                    d[f'inspiration_mol_{i}'] = hit.mol
                except FailedToAssignBondOrders as e:
                    mout.warning(f'Using hit.ligand_group.rdkit_mol {p.name=} {hit._pose_name=}')
                    d[f'inspiration_mol_{i}'] = hit.ligand_group.rdkit_mol

            data.append(d)

        return pd.DataFrame(data)

    def prepare_for_scoring(self, hit_features):

        self_features = self.products.get_present_features()

        num_in_hits_but_not_self = len(hit_features - self_features)
        self._hit_feature_coverage = (len(hit_features)-num_in_hits_but_not_self)/len(hit_features)

        self._num_new_features = len(self_features - hit_features)

        self._avg_fragmenstein_ddG = np.mean([p._fragmenstein_ddG for p in self.products])
        self._std_fragmenstein_ddG = np.std([p._fragmenstein_ddG for p in self.products])
            
        p_with_RMSD = [p._fragmenstein_mRMSD for p in self.products if p._fragmenstein_mRMSD is not None]
        if p_with_RMSD:
            self._avg_fragmenstein_mRMSD = np.mean(p_with_RMSD)
            self._std_fragmenstein_mRMSD = np.std(p_with_RMSD)
        else:
            self._avg_fragmenstein_mRMSD = None
            self._std_fragmenstein_mRMSD = None
