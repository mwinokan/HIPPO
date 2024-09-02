import mout

from collections.abc import MutableSet

from .tools import df_row_to_dict
from .compound import Compound
from .pose import Pose
from .tset import TagSet
from .pset import PoseSet

import random
import molparse as mp
import pandas as pd


class CompoundSet(MutableSet):

    ### DUNDERS

    def __init__(self, name, compounds=(), immutable=False, duplicate="warn"):

        assert isinstance(name, str)

        self._name = name
        self._elements = []
        self._immutable = immutable

        for comp in compounds:
            self.add(comp, duplicate=duplicate)

        self._fingerprint_df = None
        self._metadata_cols = []

    def __repr__(self):
        return f'CompoundSet("{self.name}", #compounds={self.num_compounds}, #poses={self.num_poses})'

    def __iter__(self):
        return iter(self.compounds)

    def __contains__(self, compound):
        if isinstance(compound, str):
            return compound in [c.name for c in self.compounds]
        else:
            return compound in self.compounds

    def __len__(self):
        return len(self.compounds)

    def __getitem__(self, key):

        if isinstance(key, slice):
            return self._elements[key]

        if isinstance(key, list):
            return CompoundSet(
                "queried", [c for c in [self[k] for k in key] if c is not None]
            )

        if isinstance(key, str):
            return [comp for comp in self if comp.name == key][0]

        # if isinstance(key,Compound):
        #     # key = str(key)
        #     matches = [comp for comp in self if comp == key][0]

        elif isinstance(key, int):
            return self._elements[key]
            # return [c for c in self][key]

        matches = [comp for comp in self if comp == key]
        # matches = [comp for comp in self if str(comp) == key]

        if len(matches) < 1:
            mout.error(f"{key} not in {self}")
            return None
        elif len(matches) > 1:
            mout.error(f"Multiple {key} in {self}")
            exit()
            return None

        return matches[0]

    def __setitem__(self, key, value):
        assert isinstance(key, int)
        [c for c in self][key] = value

    def __delitem__(self, key):
        return self.discard(key)

    ### FACTORIES

    @classmethod
    def from_df(cls, name, df, protein, verbosity=1):

        self = cls.__new__(cls)

        id_col = df.columns[0]
        mol_col = df.columns[1]

        compounds = []
        for index, row in df.iterrows():
            mol_name = row[id_col]
            if verbosity:
                mout.progress(
                    index,
                    len(df),
                    prepend="mol --> Compound",
                    append=mol_name,
                    fill="🦛",
                )
            mol = row[mol_col]
            compound = Compound.from_rdkit_mol(mol_name, mol)
            compound._set_name = name
            compound._protein_system = protein
            compounds.append(compound)

        if verbosity:
            mout.finish()

        self.__init__(name, compounds)

        return self

    @classmethod
    def from_bound_pdbs(
        cls, name, pdbs, metadata_df, pdb_pattern, tags=None, animal=None
    ):

        import os

        self = cls.__new__(cls)

        self.__init__(name)

        tags = tags or []
        tags = TagSet(tags)

        for pdb in pdbs:
            pose_name = os.path.basename(pdb).removesuffix("_bound.pdb")

            metadata = df_row_to_dict(
                metadata_df[metadata_df["crystal_name"] == pose_name]
            )

            comp_name = metadata["RealCrystalName"][-5:]

            if comp_name not in self:

                compound = Compound.from_bound_pdb(
                    comp_name, pdb, metadata, tags=tags, animal=animal
                )
                compound._set_name = name
                self.add(compound)

            else:

                compound = self[comp_name]
                pose = Pose.from_bound_pdb(compound, pdb, metadata, tags=tags)
                compound.add_pose(pose)

            ### add a pose

        return self

    ### PROPERTIES

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, a):
        self._name = a

    @property
    def compounds(self):
        return self._elements

    @property
    def num_compounds(self):
        return len(self.compounds)

    @property
    def fingerprinted_poses(self):
        return PoseSet([p for p in self.poses if p._fingerprint is not None])

    @property
    def fingerprints(self):
        return [p.fingerprint for p in self.poses if p._fingerprint is not None]

    @property
    def fingerprint_df(self):
        if self._fingerprint_df is None:
            import pandas as pd

            fingerprints = self.fingerprints
            if len(fingerprints) < 1:
                mout.error(f"no fingerprints for {self}")
            self._fingerprint_df = pd.DataFrame(fingerprints)
        return self._fingerprint_df

    @property
    def df(self):
        data = []

        for comp in self:

            d = {}

            d["_Name"] = comp.name
            d["ROMol"] = comp.mol

            d.update(comp.dict)

            data.append(d)

        return pd.DataFrame(data)

    @property
    def immutable(self):
        return self._immutable

    @immutable.setter
    def immutable(self, b):
        self._immutable = b

    @property
    def poses(self):
        return PoseSet(sum([c.poses for c in self], PoseSet()))

    @property
    def num_poses(self):
        return sum([c.num_poses for c in self])

    @property
    def reactions(self):
        return sum([[r for r in c.reactions] for c in self], [])

    @property
    def smiles(self):
        return [c.smiles for c in self]

    @property
    def names(self):
        return [c.name for c in self]

    # @property
    # def building_blocks(self):
    # 	return self.get_building_blocks()

    @property
    def dict(self):
        d = dict(
            name=self.name,
            num_compounds=len(self),
            # total_bb_amount=self.total_bb_amount,
            price=self.get_price(),
            # max_lead_time=self.max_lead_time,
            # avg_lead_time=self.avg_lead_time,
            # hit_feature_coverage=self.hit_feature_coverage,
            # num_new_features=self.num_new_features,
            # avg_fragmenstein_ddG=self.avg_fragmenstein_ddG,
            # std_fragmenstein_ddG=self.std_fragmenstein_ddG,
            # avg_fragmenstein_mRMSD=self.avg_fragmenstein_mRMSD,
            # std_fragmenstein_mRMSD=self.std_fragmenstein_mRMSD,
            # lists
            compounds=[c.dict for c in self],
            # products = [p.dict for p in self.products],
        )

        try:
            d["num_products"] = len(self.products)
            d["products"] = [p.dict for p in self.products]

        except AttributeError as e:
            mout.error(e)

        return d

    ### METHODS

    def pop(self):
        assert not self.immutable
        return self._elements.pop()
        # last = self._elements[-1]
        # self._elements = self._elements[:-1]
        # return last

    def discard(self, key):
        assert not self.immutable
        if key in self:
            i = self._elements.index(key)
            del self._elements[i]
        else:
            raise ValueError(f"{key} not in {self}")

    def remove(self, key):
        assert not self.immutable
        if key in self:
            i = self._elements.index(key)
            del self._elements[i]
        else:
            raise ValueError(f"{key} not in {self}")

    def add(self, compound, duplicate="error"):
        assert not self.immutable

        if duplicate == "no_check":
            self._elements.append(compound)
            compound._set_name = self.name

        elif compound not in self._elements:
            self._elements.append(compound)
            compound._set_name = self.name

        elif duplicate == "error":
            mout.error(f"Compound already in {self}")
            mout.var("compound", f"{compound}")
            mout.var("match", f"{self[compound]}")
            raise ValueError(f"{compound} already in {self}")

        elif duplicate == "skip":
            mout.warning(f"Skipping {compound} already in {self}")

        elif duplicate == "quiet":
            return

    def shuffle(self):
        random.shuffle(self._elements)

    # def get_random(size=1):
    #     return random.sample(self._elements,size)

    def get_present_features(self):

        features = set()
        for col in [
            c for c in self.fingerprint_df.columns if c not in self._metadata_cols
        ]:
            covered = any(self.fingerprint_df[col].values)
            if covered:
                features.add(col)

        return features

    def get_price(self):
        price = 0
        for comp in self:
            price += comp.get_price(comp.amount)
        return price

    def write(self, directory, extra_data=None, overwrite=False):

        if not overwrite:
            assert not (directory / f"{self.name}.pickle").exists(), list(
                directory.glob(f"{self.name}.pickle")
            )

        mp.write(directory / f"{self.name}.pickle", self)

        json_data = self.dict

        if extra_data:
            json_data.update(extra_data)

        mp.write(directory / f"{self.name}.json", json_data)

    def summary(self):

        mout.header(self)

        print_data = []
        for comp in self:
            print_data.append(
                dict(
                    smiles=comp.smiles,
                    set_name=comp.set_name,
                    is_pains=comp.is_pains,
                    building_blocks=comp.building_blocks,
                    cost_range_str=comp.cost_range_str,
                    lead_time=comp.lead_time,
                )
            )

        print(pd.DataFrame(print_data))

    # def get_building_blocks(self,purchaseable_only=False):

    #     # from .bb import BuildingBlockSet

    #     # bb_set = BuildingBlockSet()

    #     bb_set = CompoundSet(f'{self.name} BBs')

    #     for comp in self:
    #         if comp.building_blocks is not None:
    #             for bb in comp.building_blocks:
    #                 if not purchaseable_only or bb.purchaseable:

    #                     if bb in bb_set:

    #                         bb_set[bb].amount += 1

    #                     else:

    #                         bb = bb.copy()
    #                         bb.amount = 1
    #                         bb_set.add(bb, error_duplicate=False)

    #     return bb_set

    def get_by_smiles(self, smiles):
        matches = [c for c in self if c.smiles == smiles]
        assert len(matches) == 1
        return matches[0]

    def get_by_name(self, name):
        matches = [c for c in self if c.name == name]
        assert len(matches) == 1
        return matches[0]

    ### PROTECTION

    def __le__(self):
        raise NotImplementedError("CompoundSet.__le__")

    def __lt__(self):
        raise NotImplementedError("CompoundSet.__lt__")

    def __eq__(self):
        raise NotImplementedError("CompoundSet.__eq__")

    def __ne__(self):
        raise NotImplementedError("CompoundSet.__ne__")

    def __gt__(self):
        raise NotImplementedError("CompoundSet.__gt__")

    def __ge__(self):
        raise NotImplementedError("CompoundSet.__ge__")

    def __and__(self):
        raise NotImplementedError("CompoundSet.__and__")

    def __or__(self):
        raise NotImplementedError("CompoundSet.__or__")

    # def __sub__(self):
    #     raise NotImplementedError('CompoundSet.__sub__')

    def __xor__(self):
        raise NotImplementedError("CompoundSet.__xor__")

    def isdisjoint(self):
        raise NotImplementedError("CompoundSet.isdisjoint")

    def __ior__(self):
        raise NotImplementedError("CompoundSet.__ior__")

    def __ixor__(self):
        raise NotImplementedError("CompoundSet.__ixor__")

    def __iand__(self):
        raise NotImplementedError("CompoundSet.__iand__")

    def __isub__(self):
        raise NotImplementedError("CompoundSet.__isub__")
