"""
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
"""

import os

if __name__ == '__main__':
	base = os.path.dirname(os.path.realpath(__file__))
	files = os.listdir(base)
	tests = [f for f in files if '.py' in f and 'swigtest.py' != f and 'test_' in f]
	for test in tests:
		print(14*'*****'+'\nRunning: '+test+'\n'+14*'*****'+'\n')
		os.system('python '+base+'/'+test)
