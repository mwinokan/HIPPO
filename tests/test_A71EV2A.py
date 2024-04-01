#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.path.abspath("../"))

import unittest
from hippo import HIPPO

class TestA71EV2A(unittest.TestCase):

	def setUp(self):
		os.system('tar -xf ../data/A71EV2A.tar.gz')
		os.system('rm test_A71EV2A.sqlite')
		self.animal = HIPPO('TestA71EV2A', 'test_A71EV2A.sqlite')

	def test_add_hits(self):
		self.animal.add_hits(
			target_name='A71EV2A',
			metadata_csv='A71EV2A/metadata.csv',
			aligned_directory='A71EV2A/aligned_files',
		)

		self.assertEqual(len(self.animal.compounds),84)
		self.assertEqual(len(self.animal.poses),107)

	def tearDown(self):
		self.animal.db.close()

if __name__ == '__main__':
	unittest.main()
