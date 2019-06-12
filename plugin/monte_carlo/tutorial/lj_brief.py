import feasst
feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.add(feasst.Configuration(feasst.args(
  {"cubic_box_length": "8", "particle_type": "../../../forcefield/data.lj"})))
mc.add(feasst.Potential(feasst.MakeModelLJ()))
mc.add(feasst.Potential(feasst.MakeLongRangeCorrections()))
mc.add(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": "1.5", "chemical_potential": "1."})))
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "tunable_param": "2."})))
mc.add(feasst.MakeLog(feasst.args({"steps_per" : "10000"})))
mc.add(feasst.MakeMovie(feasst.args(
  {"steps_per" : "10000", "file_name" : "movie.xyz"})))
mc.add(feasst.MakeCheckEnergy(feasst.args(
  {"steps_per" : "10000", "tolerance" : "1e-10"})))
mc.add(feasst.MakeTuner(feasst.args({"steps_per" : "10000"})))
mc.seek_num_particles(50)
mc.attempt(int(1e6))
