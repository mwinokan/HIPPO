#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO


class TestCHIKV_Mac(unittest.TestCase):

    def setUp(self):
        # os.system("tar -xf ../data/CHIKV_Mac.tar.gz")

        from hippo.fragalysis import download_target

        download_target("CHIKV_Mac", destination="targets", overwrite=True, unzip=True)

        os.system("rm test_CHIKV_Mac.sqlite")
        self.animal = HIPPO("TestCHIKV_Mac", "test_CHIKV_Mac.sqlite")

    def test_add_hits(self):
        self.animal.add_hits(
            target_name="CHIKV_Mac",
            metadata_csv="targets/CHIKV_Mac/metadata.csv",
            aligned_directory="targets/CHIKV_Mac/aligned_files",
        )

        # self.assertEqual(len(self.animal.compounds), 84)
        # self.assertEqual(len(self.animal.poses), 107)

    def tearDown(self):
        self.animal.db.close()


if __name__ == "__main__":
    unittest.main()
