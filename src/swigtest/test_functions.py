"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
"""

"""@docstring
@brief Tests functions
@author Nathan A. Mahynski
@date 5/16/2017
@filename test_functions.py
"""

import sys, os
import unittest
import numpy as np

class TestFunctions(unittest.TestCase):
  #check that random number generator selects particles randomly
  def testMyRanNum(self):
    imin=5
    imax=9
    n=10000
    x=np.array(range(imax+1))
    feasst.ranInitByDate()
    s=feasst.Space(3,0)
    for i in range(n):
        ran=s.uniformRanNum(imin, imax)
        x[ran]=x[ran]+1
    for i in range(imin,imax+1):
        self.assertAlmostEqual(1./(imax-imin+1),x[i]/float(n), 1)

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
