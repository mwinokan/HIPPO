
from .set import CompoundSet
import pandas as pd
import mout
import json
import numpy as np
from collections.abc import Set
from pprint import pprint

class BuildingBlock:

    ### DUNDERS

    def __init__(self, smiles, molecular_weight=None, catalog_metadata=None, minimise='cost_max'):
        
        self._smiles = smiles
        self._molecular_weight = molecular_weight
        self._cost_str = None
        self._cost_min = None
        self._cost_max = None
        self._cost_unit = None
        self._lead_time = None
        self._catalog_metadata = None

        self.sanitise_catalog_metadata(catalog_metadata)

    def __repr__(self):
        if self._cost_str is not None:
            return f'BuildingBlock({self.smiles}, cost={self._cost_str}, lead_time={self.lead_time} weeks)'
        else:
            return f'BuildingBlock({self.smiles})'

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __hash__(self):
        return hash(self.smiles)

    ### PROPERTIES

    @property
    def smiles(self):
        return self._smiles

    @property
    def lead_time(self):
        return self._lead_time

    ### METHODS

    def get_purchaseable_entries(self):
        return [ data for data in self._catalog_metadata if 'purchaseInfo' in data]

    # def get_exact_catalogue_entries(self):
    #     return [ data for data in self.get_purchaseable_entries() if data['inchikeyMatches']['exact']]

    def get_shortest_lead_time(self):
        
        entries = self.get_purchaseable_entries()

        if not entries:
            mout.error(f'No purchaseable entries in catalog_metadata ({self})')
            return

        self._best_catalog_entry = sorted(entries, key=lambda x: x['purchaseInfo']['lead_time_str'])[0]

        self._cost_str = self._best_catalog_entry['purchaseInfo']['price_range_str']
        self._cost_min = self._best_catalog_entry['purchaseInfo']['price_range_min']
        self._cost_max = self._best_catalog_entry['purchaseInfo']['price_range_max']
        self._cost_unit = self._best_catalog_entry['purchaseInfo']['price_range_unit']
        self._lead_time = self._best_catalog_entry['purchaseInfo']['lead_time_str']

    def get_cost(self, field='max'):
        
        if self._cost_min == 0:
            return self._cost_max

        if self.cost_max == np.inf:
            return self._cost_min

        if field == 'min':
            return self._cost_min

        if field == 'max':
            return self._cost_max

    def get_lowest_cost(self, use_max=False):
        
        entries = self.get_purchaseable_entries()

        if not entries:
            mout.error(f'No purchaseable entries in catalog_metadata ({self})')
            return

        if use_max:
            key = 'price_range_max'
        else:
            key = 'price_range_min'

        self._best_catalog_entry = sorted(entries, key=lambda x: x['purchaseInfo'][key])[0]

        self._cost_str = self._best_catalog_entry['purchaseInfo']['price_range_str']
        self._cost_min = self._best_catalog_entry['purchaseInfo']['price_range_min']
        self._cost_max = self._best_catalog_entry['purchaseInfo']['price_range_max']
        self._cost_unit = self._best_catalog_entry['purchaseInfo']['price_range_unit']
        self._lead_time = self._best_catalog_entry['purchaseInfo']['lead_time_str']

    def sanitise_catalog_metadata(self, catalog_metadata):

        if catalog_metadata is None:
            return

        try:
            self._catalog_metadata = eval(catalog_metadata)
        except TypeError as e:
            mout.error(f'{self}: {e}')
            pprint(catalog_metadata)
            exit()
        # self._catalog_metadata = json.loads(catalog_metadata.strip('"').replace("'",'"'))

        for entry in self.get_purchaseable_entries():

            if entry['purchaseInfo']['isBuildingBlock']:
                price_range_str = entry['purchaseInfo']['bbPriceRange']
                lead_time_str = entry['purchaseInfo']['bbLeadTimeWeeks']
            else:
                price_range_str = entry['purchaseInfo']['scrPriceRange']
                lead_time_str = entry['purchaseInfo']['scrLeadTimeWeeks']

            if price_range_str == 'unknown':
                del entry['purchaseInfo']
                continue

            price_range_str = price_range_str.replace('k','000')
            
            split_price_range = price_range_str.split('/')
            price_range_unit = f'/{split_price_range[1].strip()}'

            if '-' in price_range_str:
                price_range_min = float(split_price_range[0].split('-')[0].lstrip('$').strip())
                price_range_max = float(split_price_range[0].split('-')[1].strip())
            elif '>' in price_range_str:
                price_range_max = np.inf
                price_range_min = float(split_price_range[0].lstrip('> $').strip())
            elif '<' in price_range_str:
                price_range_min = 0
                price_range_max = float(split_price_range[0].lstrip('< $').strip())
            else:
                mout.error(price_range_str)
                raise Exception('Could not parse price_range_str')

            # if price_range_unit == '/mg':
            #     price_range_min *= 1000
            #     price_range_max *= 1000
            #     price_range_unit = '/g'

            entry['purchaseInfo']['price_range_str'] = price_range_str
            entry['purchaseInfo']['price_range_min'] = price_range_min
            entry['purchaseInfo']['price_range_max'] = price_range_max
            entry['purchaseInfo']['price_range_unit'] = price_range_unit
            entry['purchaseInfo']['lead_time_str'] = lead_time_str

class BuildingBlockSet(Set):

    ### DUNDERS

    def __init__(self, bbs=()):
        # super(BuildingBlockSet,self).__init__(s)
        self._name = None
        self._cost = None

        self._elements = []

        for bb in bbs:
            if bb not in self._elements:
                self._elements.append(bb)

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
            return f'BuildingBlockSet(#building_blocks={len(self)})'

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def cost(self):
        return self._cost

    @property
    def smiles(self):
        return[ bb.smiles for bb in self._elements]
    
    ### METHODS

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

        return comp_set

    def summary(self):

        mout.header(self)

        print_data = []
        for bb in self:
            print_data.append(dict(
                smiles=bb.smiles,
                cost=bb._cost_str,
                lead_time=bb.lead_time,
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
