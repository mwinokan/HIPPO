#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO


class TestCHIKV_Mac(unittest.TestCase):

    def setUp(self):
        self.animal = HIPPO("TestCHIKV_Mac", "test_CHIKV_Mac.sqlite")

        self.animal.db.insert_subsite(target=1, name="ribose")
        self.animal.db.insert_subsite(target=1, name="adenine")

    def test_web(self):

        from hippo.web import ProposalPage

        page = ProposalPage(
            output_dir="/Users/tfb64483/Software/HIPPO/tests",
            animal=self.animal,
        )

        page.write_html()

    def tearDown(self):
        self.animal.db.close()


if __name__ == "__main__":
    unittest.main()
