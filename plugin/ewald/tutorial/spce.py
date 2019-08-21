import feasst

def lrc():
  potential = feasst.Potential()
  potential.set_visit_model(feasst.MakeLongRangeCorrections())
  return potential

def system(config, alphaL = 5.6, kmax_squared=27):
  minl = config.domain().min_side_length()
  config.add_model_param("alpha", alphaL/minl)
  system = feasst.System()
  feasst.add_ewald_with(feasst.MakeModelLJ(), config, system, kmax_squared)
  system.add(lrc())
  system.add(config)
  system.precompute()
  return system
