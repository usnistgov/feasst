#! /usr/bin/env python
import os, sys
feasstdir = os.getenv("HOME") + "/feasst"
sys.path.append(feasstdir + "/src")
import feasst, pyfeasst
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

