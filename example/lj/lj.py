#! /usr/bin/env python
import os, sys
feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/build"
if (not os.path.isfile(feasstdir+"/_feasst.so")):
  feasstdir = os.getenv("FEASST_INSTALL_DIR_") + "/src"
sys.path.append(feasstdir)
import feasst
space = feasst.Space(3, 0)
for dim in range(space.dimen()): space.lset(8, dim) # 8 box length
space.addMolInit("../../forcefield/data.lj")
pair = feasst.PairLJ(space, 3)    # potential truncation at 3
pair.initEnergy()
criteria = feasst.CriteriaMetropolis(1.2, 1.);  # 1/kT = 1.2
mc = feasst.MC(space, pair, criteria)
maxMoveParam = 0.1
feasst.transformTrial(mc, "translate", maxMoveParam)
mc.nMolSeek(50, "../../forcefield/data.lj")   # add 50 particles
mc.initLog("log", int(1e4))
mc.initMovie("movie", int(1e4))
mc.runNumTrials(int(1e6))
# mc.initWindows(1)       # initialize parallel windows
# mc.runNumSweeps(1, -1)  # run until all windows sweep once

