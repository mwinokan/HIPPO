#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO

os.system("rm -v test_radical.sqlite")
animal = HIPPO("TestRadical", "test_radical.sqlite")

from mlog import setup_logger

logger = setup_logger("TestRadical")


class TestRadical(unittest.TestCase):

    def test_merge_ok(self):
        logger.header("test_merge_ok")
        compound = animal.register_compound(
            smiles="Nc1ccc(S(=O)(=O)Nc2cc3c(OCCO3)c(CO)n2)cc1", radical="error"
        )
        print(compound)
        assert compound

    def test_merge_bad_aromaticity(self):
        logger.header("test_merge_bad_aromaticity")
        compound = animal.register_compound(
            smiles="Nc1ccc(S(=O)(=O)NC2CC3=C(OCCO3)C(CO)N2)cc1", radical="error"
        )
        print(compound)
        assert compound

    def test_merge_bad_aromaticity_and_radicals(self):
        logger.header("test_merge_bad_aromaticity_and_radicals")
        compound = animal.register_compound(
            smiles="Nc1ccc(S(=O)(=O)NC2CC3=C(O[CH][CH]O3)C(CO)N2)cc1", radical="error"
        )
        print(compound)
        assert not compound

    def test_merge_fix_radicals(self):
        logger.header("test_merge_fix_radicals")
        compound = animal.register_compound(
            smiles="Nc1ccc(S(=O)(=O)NC2CC3=C(O[CH][CH]O3)C(CO)N2)cc1", radical="remove"
        )
        print(compound)
        assert compound


if __name__ == "__main__":
    unittest.main()
    animal.db.close()
