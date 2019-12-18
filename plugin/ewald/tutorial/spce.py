import feasst

def system(config, alphaL = 5.6, kmax_squared=38, rcut=False):
    config.add_model_param("alpha", alphaL/config.domain().min_side_length())
    if rcut:
        config.set_model_param("cutoff", 0, rcut)
        config.set_model_param("cutoff", 1, rcut)
    system = feasst.System()
    system.add(config)
    feasst.add_ewald_with(feasst.MakeLennardJones(), system, kmax_squared)
    system.add(feasst.Potential(feasst.MakeLongRangeCorrections()))
    # system.precompute()
    return system
