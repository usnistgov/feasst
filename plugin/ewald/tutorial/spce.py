import feasst

def lrc():
    potential = feasst.Potential()
    potential.set_visit_model(feasst.MakeLongRangeCorrections())
    return potential

def system(config, alphaL = 5.6, kmax_squared=27):
    config.add_model_param("alpha", alphaL/config.domain().min_side_length())
    system = feasst.System()
    system.add(config)
    feasst.add_ewald_with(feasst.MakeModelLJ(), config, system, kmax_squared)
    #system.add(lrc())
    #system.precompute()
    return system
