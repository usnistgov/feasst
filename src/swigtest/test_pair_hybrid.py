"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
"""

"""@docstring
@brief Tests pair interactions
@author Nathan A. Mahynski and Harold Hatch
@date 5/16/2017
@filename test_pair_hybrid.py
"""

import sys, os
import unittest
import numpy as np

class TestPairRestart(unittest.TestCase):
	def testLJMultiRcut(self):
		s = feasst.Space(3,0)
		s.lset(10)
		s.addMolInit('../unittest/data.lj')
		s.addMolInit('../unittest/data.ljb')

		ipair = feasst.PairHybrid (s, '../unittest/rstpair')
		ipair.initEnergy()

		self.assertTrue (ipair.getPair(0).rCutij(0,0) == (ipair.getPair(0).sig(0)+ipair.getPair(0).sig(0))/2*3.0)
		self.assertTrue (ipair.getPair(0).rCutij(0,1) == (ipair.getPair(0).sig(0)+ipair.getPair(0).sig(1))/2*3.0)
		self.assertTrue (ipair.getPair(0).rCutij(1,1) == (ipair.getPair(0).sig(1)+ipair.getPair(0).sig(1))/2*3.0)

if __name__ == "__main__":
	def find_all(name, path):
	    result = []
	    for root, dirs, files in os.walk(path):
	        if name in files:
	            result.append(os.path.join(root, name))
	    return result

	feasstname = find_all("_feasst.so", os.path.dirname(os.path.realpath(__file__))+"/../../")
	if (len(feasstname) > 1):
		print ('Found multiple FEASST sources which may cause an error: ', feasstname)
		sys.exit()
	else:
		feasstdir = feasstname[0].split('/_feasst.so')[0]
		sys.path.append(feasstdir)

	import feasst
	unittest.main()
