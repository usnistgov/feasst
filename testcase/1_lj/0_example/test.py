import feasst
space = feasst.makeSpace(feasst.args(
    {"dimen" : "3",                 # 3D space
     "boxLength" : "8"}))           # cubic periodic boundaries
pair = feasst.makePairLJ(space, feasst.args(
    {"rCut" : "3",                  # potential truncation
     "cutType" : "lrc"}))           # long range corrections
criteria = feasst.makeCriteriaMetropolis(feasst.args(
    {"beta" : "1.2"}))              # beta = 1/k_B/T
mc = feasst.MC(pair, criteria)
feasst.addTrialTransform(mc, feasst.args(
    {"type" : "translate",          # attempt particle translations
     "maxMoveParam" : str(0.1)}))   # maximum displacement for each dimension
mc.nMolSeek(50)                     # add particles
mc.initLog("log", int(1e4))         # output instantaneous values
mc.initMovie("movie", int(1e4))     # output xyz trajectory
mc.runNumTrials(int(1e6))           # perform MC trials
