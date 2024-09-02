#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO

os.system("rm -v test_recipe1.sqlite")
animal = HIPPO("TestRecipe1", "test_recipe1.sqlite")


class TestRecipe1(unittest.TestCase):

    def register_compounds(self):

        smiles = [
            "O=C(O)C1CCCO1",
            "Cn1ncc(NC(=O)OC(C)(C)C)c1N",
            "N#CCC(=O)O",
            "CC1(C(=O)O)CCCO1",
            "O=C1CCC(C(=O)O)O1",
            "CC1CCC(C(=O)O)O1",
            "O=C(O)C1OCCC1F",
            "CCn1ncc(NC(=O)OC(C)(C)C)c1N",
            "CC1COC(C(=O)O)C1",
            "O=C1COC(C(=O)O)C1",
        ]

        for s in smiles:
            animal.register_compound(smiles=s)

        assert len(animal.compounds) == 10

    def register_reactions(self):

        animal.register_reaction(type="combination", product=1, reactants=[2, 3])
        animal.register_reaction(type="deprotection", product=3, reactants=[4])
        animal.register_reaction(type="deprotection", product=3, reactants=[7])
        animal.register_reaction(type="combination", product=4, reactants=[5, 6])
        animal.register_reaction(type="combination", product=7, reactants=[8, 9])
        animal.register_reaction(type="deprotection", product=7, reactants=[8])

        assert len(animal.reactions) == 6

    def test_recipe(self):

        self.register_compounds()
        self.register_reactions()

        recipes = animal.reactions[0].get_recipe(debug=False, pick_cheapest=False)

        assert len(recipes) == 3

        for recipe in recipes:
            assert recipe.products[0].amount == 1
            assert len(recipe.intermediates) == 2
            assert len(recipe.reactants) <= 3
            assert len(recipe.products) == 1


if __name__ == "__main__":
    unittest.main()
    animal.db.close()
