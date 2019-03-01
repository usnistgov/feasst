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

def translate():
  trial = feasst.TrialTranslateShrPtr()
  trial.set_weight(1)
  trial.set_max_move(2.)
  return trial

mc = feasst.MonteCarlo()
mc.set(system())
mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": "1.2", "add_activity": "1."})))
mc.add(translate())
mc.seek_num_particles(50);
num_periodic = int(1e4);
mc.add(feasst.MakeLog(feasst.args(
  {"steps_per" : str(num_periodic)})))
mc.add(feasst.MakeMovie(feasst.args(
  {"steps_per" : str(num_periodic),
   "file_name" : "tmp/lj50movie.xyz"})));
mc.add(feasst.MakeEnergyCheck(feasst.args(
  {"steps_per" : str(num_periodic),
   "tolerance" : "1e-10"})));
mc.add(feasst.MakeTuner(feasst.args(
  {"steps_per" : str(num_periodic)})));
mc.attempt(int(1e6))

print("hi", "box length", mc.system().configuration().domain().side_length().str())
