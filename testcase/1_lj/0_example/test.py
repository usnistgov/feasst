import feasst
feasst.ranInitByDate()
space = feasst.makeSpace(3)
space.initBoxLength(8)
space.addMolInit(space.install_dir() + "/forcefield/data.lj")
pair = feasst.PairLJ(space, 3, # potential truncation at 3
    feasst.args({"cutType" : "lrc"}))
criteria = feasst.CriteriaMetropolis(1.2,  # beta = 1/k_B/T
    1.)  # lnz = beta mu - 3ln(Lambda)
mc = feasst.MC(space, pair, criteria)
maxMoveParam = 0.1  # maximum displacement for each dimension
feasst.transformTrial(mc, "translate", maxMoveParam)
mc.nMolSeek(50)   # add particles
mc.initLog("log", int(1e4))
mc.initMovie("movie", int(1e4))
mc.runNumTrials(int(1e6))
