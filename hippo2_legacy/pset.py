# set of Poses

import pandas as pd
from collections.abc import MutableSet


class PoseSet(MutableSet):

    def __init__(self, poses=(), immutable=False):

        self._elements = []
        self._immutable = immutable
        self._fingerprint_df = None
        self._metadata_cols = []

        for pose in poses:
            self.add(pose)

    ### FACTORIES

    ### PROPERTIES

    @property
    def immutable(self):
        return self._immutable

    @immutable.setter
    def immutable(self, b):
        self._immutable = b

    @property
    def poses(self):
        return self._elements

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
    def df(self, compound_dict=True):
        data = []

        for pose in self:

            d = {}

            d["_Name"] = pose.longname
            d["ROMol"] = pose.mol

            if compound_dict:
                d["compound"] = pose.compound.dict
                d["compound_mol"] = pose.compound.mol

            if compound_dict:
                d["base"] = pose.compound.base.dict
                d["base_mol"] = pose.compound.base.mol

            d.update(pose.dict)

            data.append(d)

        return pd.DataFrame(data)

    ### METHODS

    def pop(self):
        assert not self.immutable
        return self._elements.pop()

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

    def add(self, compound):
        assert not self.immutable
        # if compound not in self._elements:
        self._elements.append(compound)
        # else:
        # raise ValueError(f'{compound} already in {self}')

    def remove_unfingerprinted(self, animal):

        new = []
        for pose in self:
            if pose.fingerprint is not None:
                new.append(pose)

        self.__init__(new)

    def get_present_features(self):

        features = set()
        for col in [
            c for c in self.fingerprint_df.columns if c not in self._metadata_cols
        ]:
            covered = any(self.fingerprint_df[col].values)
            if covered:
                features.add(col)

        return features

    ### DUNDERS

    def __contains__(self, pose):
        if isinstance(pose, str):
            return pose in [c.name for c in self.poses]
        else:
            return pose in self.poses

    def __repr__(self):
        if not self:
            return f"PoseSet(empty)"
        else:
            return f'PoseSet(#poses={len(self)}, [{", ".join(p.name for p in self)}])'

    def __len__(self):
        return len(self._elements)

    def __iter__(self):
        return iter(self._elements)

    def __add__(self, other):
        if isinstance(other, PoseSet):
            return PoseSet(self._elements + other._elements)
        elif isinstance(other, list):
            return PoseSet(self._elements + other)

    def __iadd__(self, other):
        if isinstance(other, PoseSet):
            return PoseSet(self._elements + other._elements)
        elif isinstance(other, list):
            return PoseSet(self._elements + other)

    def __getitem__(self, key):
        return self._elements[key]

    def __iter__(self):
        return iter(self.poses)
