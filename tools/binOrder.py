#! /usr/bin/env python

from libxdrfile import xdrfile_open, xdrfile_close, read_xtc_natoms, read_xtc, DIM, exdrOK
import numpy as np
import math, os, sys
feasstdir = os.getenv("HOME") + "/feasst"
sys.path.append(feasstdir + "/src")
sys.path
import feasst

# initialize cell list and cluster definition
s = feasst.Space("tmp/rstspace")
s.initCellAtomCut(1)
# obtain rCut for cluster definitions from pair restart file
pfactory = feasst.PairLJ(s, 0)  # temporary generic pair class
pair = pfactory.makePair(s, "tmp/rstpair")
s.updateCells(pair.rCut())
s.addTypeForCluster(1)
s.preMicellarAgg(5)

# initialize histogram for average cluster size as a function of order parameter
clusterAccVec = feasst.AccumulatorVec()

#for p in range(1, 3):
for p in range(0, 12):

  # loop through xtc file and log file
  # obtain the order parameter from the log file
  f = open('logp'+str(p)+'pr', 'r')
  lines = f.readlines()
  print len(lines)
  # open xtc file
  xtcFileName = 'moviep'+str(p)+'prn75.xtc'
  xtcFile = xdrfile_open(xtcFileName, 'r')
  endXTC = s.readXTC(xtcFileName, xtcFile)
  conf = 0
  while ( (endXTC == 0) and (conf != len(lines)) ):
   
    # compute clusters
    s.wrapMol()
    s.updateCellofallMol()
    s.updateClusters(pair.rCut())
    s.xClusterGen()
    #s.xClusterShape()
   
    # obtain order param
    columns = lines[conf].split(' ')
    currentOrder = float(columns[3])
  
    # store stats and prep for next conf
    clusterAccVec.accumulate(s.nClusters(), currentOrder)
    #print "nMol", s.nMol(), "nClusters", s.nClusters(), "order", currentOrder
    endXTC = s.readXTC(xtcFileName, xtcFile)
    conf = conf + 1

# obtain final stats
plotx = list()
ploty = list()
for i in range(clusterAccVec.size()):
  if (clusterAccVec.vec(i).nValues() > 20):
    plotx.append(i)
    ploty.append(clusterAccVec.vec(i).average())

# plot
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
plt.plot(plotx, ploty)
plt.xlabel("nCluster")
plt.ylabel("average order")
plt.savefig(os.path.splitext(__file__)[0] + ".eps")

# finally close file
s.cellOff() # avoid error check
xdrfile_close(xtcFile)
