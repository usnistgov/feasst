"""@docstring
@brief Tests Monte Carlo
@author Nathan A. Mahynski
@date 5/16/2017
@filename test_mc.py
"""

import sys, os
import unittest
import numpy as np

class TestMC(unittest.TestCase):
  def testMC(self):
    s=feasst.Space(3,0)
    s.lset(10)
    p=feasst.PairLJMulti(s, 3)
    c=feasst.CriteriaWLTMMC(1, 0.1, "nmol",0.5,1.5,2)
    mc=feasst.WLTMMC(s, p, c)

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
