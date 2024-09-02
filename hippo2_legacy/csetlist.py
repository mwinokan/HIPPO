import mout
from collections import UserList
from .cset import CompoundSet


class CompoundSetList(UserList):

    def __init__(self, inherited_list=list()):
        super(CompoundSetList, self).__init__(inherited_list)

    def __getitem__(self, key):

        if isinstance(key, slice):
            return CompoundSetList(self.data[key])

        if isinstance(key, int):
            return self.data[key]

        if isinstance(key, str):
            for cs in self.data:
                if cs.name == key:
                    return cs
            else:
                mout.error(f"No CompoundSet named {key}")

    def __setitem__(self, key, value):

        if isinstance(key, slice):
            self.data[key] = value

        if isinstance(key, int):
            data[key] = value

        if isinstance(key, str):
            for cs in self.data:
                if cs.name == key:
                    cs = value
                    return
            else:
                mout.error(f"No CompoundSet named {key}")

    def __contains__(self, key):
        assert isinstance(key, str)
        return key in [cs.name for cs in self.data]

    @property
    def names(self):
        return [cs.name for cs in self.data]

    def append(self, item):
        if item.name in self.names:
            mout.error(
                f"CompoundSetList already contains a CompoundSet named {item.name}"
            )
            return

        self.data.append(item)

    @property
    def all_compounds(self):
        # comps = set(sum([cs.compounds for cs in self],[]))
        comps = sum([cs.compounds for cs in self], [])
        cset = CompoundSet("all_compounds", ())
        for c in comps:
            cset.add(c, duplicate="quiet")
        return cset
