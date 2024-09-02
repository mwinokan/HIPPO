# set of Reactions

from collections.abc import MutableSet


class ReactionSet(MutableSet):

    def __init__(self, reactions=(), immutable=False):

        self._elements = []
        self._immutable = immutable

        for reaction in reactions:
            self.add(reaction)

    ### FACTORIES

    ### PROPERTIES

    @property
    def reactions(self):
        return self._elements

    @property
    def immutable(self):
        return self._immutable

    @immutable.setter
    def immutable(self, b):
        self._immutable = b

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
        if compound not in self._elements:
            self._elements.append(compound)
        else:
            raise ValueError(f"{compound} already in {self}")

    def remove_unpurchaseable(self, animal):

        new = []
        for reaction in self:

            for bb in reaction.reactants:

                if not bb.price_picker:
                    break

                if bb.lead_time > animal.max_lead_time:
                    break

                if bb.get_price(animal.min_bb_quantity) > animal.max_bb_price:
                    break

            else:
                new.append(reaction)

        self.__init__(new)

    ### DUNDERS

    def __getitem__(self, key):
        return self._elements[key]
        # if isinstance(key, int):
        # raise Exception('unsupported key type')

    def __contains__(self, reactions):
        # if isinstance(reactions,str):
        # 	return reactions in [c.name for c in self.reactions]
        # else:
        return reactions in self.reactions

    def __repr__(self):
        if not self:
            return f"ReactionSet(empty)"
        else:
            return f'ReactionSet(#reactions={len(self)}, [{", ".join(r.type for r in self)}])'

    def __len__(self):
        return len(self._elements)

    def __iter__(self):
        return iter(self._elements)
