"""@docstring
@brief Tests space
@author Nathan A. Mahynski
@date 5/16/2017
@filename test_space.py
"""

import sys, os
import unittest
import numpy as np

class TestSpace(unittest.TestCase):
    #check that space class returns the correct number of atoms
    def testNatoms(self):
        for natoms in range(1,10):
            for dimen in range(1,4):
                s=feasst.Space(dimen,0)
                self.assertEqual(dimen,s.dimen())

    #check that init_config works as expected
    def testInitConfandXaccessor(self):
        for natoms in range(1,10):
            for dimen in range(1,4):
                s=feasst.Space(dimen,0)
                s.init_config(natoms)
                self.assertEqual((natoms-1)*dimen*0.95,s.x(natoms-1,0))

    #checks setting positions of particles
    def testXset(self):
        s=feasst.Space(3,0)
        s.init_config(3)
        x0=4.
        s.xset(x0,0,0)
        self.assertEqual(x0,s.x(0,0))

    #checks setting box dimensions
    def testLset(self):
        s=feasst.Space(3,0)
        l0=4.
        self.assertEqual(0.,s.l(0))
        s.lset(4.,0)
        self.assertEqual(l0,s.l(0))

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
