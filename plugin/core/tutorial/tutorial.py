import feasst

feasst.seed_random_by_date()

def config():
  config = feasst.Configuration()
  domain = feasst.Domain()
  domain.set_cubic(8.)
  config.set_domain(domain)
  config.add_particle_type("../../../forcefield/data.lj")
  return config

def lj():
  potential = feasst.Potential()
  potential.set_model(feasst.ModelLJShrPtr())
  return potential

def lrc():
  potential = feasst.Potential()
  potential.set_visit_model(feasst.LongRangeCorrectionsShrPtr())
  return potential

def system():
  system = feasst.System()
  system.add(config())
  system.add(lj())
  system.add(lrc())
  return system

def criteria():
  criteria = feasst.CriteriaMetropolis()
  criteria.set_beta(1.2)
  criteria.add_activity(1.)
  return criteria

def translate(domain):
  trial = feasst.TrialTranslateShrPtr()
  trial.set_weight(1)
  trial.set_max_move(2.)
  trial.set_max_move_bounds(domain)
  return trial

def log():
  log = feasst.LogShrPtr()
  log.set_steps_per(num_periodic)
  return log

def movie():
  movie = feasst.MovieShrPtr()
  movie.set_steps_per(num_periodic)
  movie.set_file_name("movie.xyz")
  return movie

def tune():
  tune = feasst.TunerShrPtr()
  tune.set_steps_per(num_periodic)
  return tune

def check():
  check = feasst.EnergyCheckShrPtr()
  check.set_steps_per(num_periodic)
  check.set_tolerance(1e-10)
  return check

mc = feasst.MonteCarlo()
mc.set(system())
mc.set(criteria())
mc.add(translate(mc.system().configuration().domain()))
mc.seek_num_particles(50);
num_periodic = int(1e4);
mc.add(log())
mc.add(check())
mc.add(movie())
mc.add(tune())
mc.attempt(int(1e6))

print("hi", "box length", mc.system().configuration().domain().side_length().str())
