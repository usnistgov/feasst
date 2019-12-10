import feasst

def system(config, alphaL = 5.6, kmax_squared=38):
    config.add_model_param("alpha", alphaL/config.domain().min_side_length())
    system = feasst.System()
    system.add(config)
    feasst.add_ewald_with(feasst.MakeLennardJones(), system, kmax_squared)
    system.add(feasst.Potential(feasst.MakeLongRangeCorrections()))
    # system.precompute()
    return system
