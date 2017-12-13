import feasst
space = feasst.makeSpace(3,             #3D space
    feasst.args({"boxLength" : "8"}))   # cubic box length
pair = feasst.PairLJ(space, feasst.args(
    {"rCut" : "3",          # potential truncation
     "cutType" : "lrc"}))   # long range corrections
criteria = feasst.CriteriaMetropolis(1.2)  # beta = 1/k_B/T
mc = feasst.MC(space, pair, criteria)
feasst.transformTrial(mc, feasst.args(
    {"type" : "translate",
     "maxMoveParam" : str(0.1)}))  # maximum displacement for each dimension
mc.nMolSeek(50)   # add particles
mc.initLog("log", int(1e4))
mc.initMovie("movie", int(1e4))
mc.runNumTrials(int(1e6))
