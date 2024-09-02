#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO

os.system("rm -v test_recipe2.sqlite")
animal = HIPPO("TestRecipe2", "test_recipe2.sqlite")

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

animal.register_reaction(type="deprotection", product=1, reactants=[3])
animal.register_reaction(type="combination", product=3, reactants=[4, 5])
animal.register_reaction(type="deprotection", product=3, reactants=[6])

animal.register_reaction(type="combination", product=2, reactants=[7, 8])
animal.register_reaction(type="combination", product=7, reactants=[6, 9])
animal.register_reaction(type="deprotection", product=7, reactants=[10])


class TestRecipe2(unittest.TestCase):

    def setUp(self):
        self.recipes = animal.compounds[1, 2].get_recipe(amount=3, pick_cheapest=False)
        assert len(self.recipes) == 4

    def test_register_compounds(self):
        assert len(animal.compounds) == 10

    def test_register_reactions(self):
        assert len(animal.reactions) == 6

    def test_recipe(self):

        for recipe in self.recipes:
            assert len(recipe.intermediates) == 2
            assert len(recipe.reactants) <= 5
            assert len(recipe.products) == 2

            for product in recipe.products:
                assert product.amount == 3

            for reactant in recipe.reactants:
                assert reactant.amount % 3 == 0

    def test_graph(self):

        recipe = self.recipes[0]

        fig = recipe.draw(graph_only=True)

        assert len(recipe.intermediates) == 2
        assert len(recipe.reactants) <= 5
        assert len(recipe.products) == 2

        for product in recipe.products:
            assert product.amount == 3

        for reactant in recipe.reactants:
            assert reactant.amount % 3 == 0

    def test_sankey(self):

        recipe = self.recipes[0]

        fig = recipe.sankey()

        assert len(recipe.intermediates) == 2
        assert len(recipe.reactants) <= 5
        assert len(recipe.products) == 2

        for product in recipe.products:
            assert product.amount == 3

        for reactant in recipe.reactants:
            assert reactant.amount % 3 == 0


if __name__ == "__main__":
    unittest.main()
    animal.db.close()
