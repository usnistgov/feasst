#! /usr/bin/env python

import os, sys
#from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
feasstdir = os.getenv("HOME") + "/feasst"
sys.path.append(feasstdir + "/src")
import feasst, pyfeasst
import glob, argparse
#import glob, math, argparse

parser = argparse.ArgumentParser()
parser.add_argument("--colMatName", "-i", help="collection matrix file", default="colMat", type=str)
args = parser.parse_args()

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
#fig = plt.figure(figsize=size)

#for p in range(12):
#  print p
criteria = feasst.CriteriaWLTMMC(args.colMatName)
cv = pyfeasst.DoubleVector2list(criteria.heatCapacity())
beta = []
for bin in range(criteria.nBin()):
  beta.append(criteria.bin2m(bin))
  print criteria.bin2m(bin), cv[bin]
#print beta, cv 
plt.plot(beta, cv)

plt.xlabel(r"$\beta$", fontsize=18)
plt.ylabel(r"$C_V$", fontsize=18)
plt.tight_layout()
name = os.path.splitext(os.path.basename(__file__))[0]
plt.savefig(name+".eps")
plt.savefig(name+".jpg")




