import feasst
feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.add(feasst.Configuration(feasst.args(
  {"cubic_box_length": "8", "particle_type": "../../../forcefield/data.lj"})))
mc.add(feasst.Potential(feasst.MakeModelLJ()))
mc.add(feasst.Potential(feasst.MakeLongRangeCorrections()))

#mc.add(feasst.Potential(feasst.MakeModelHardShape(feasst.MakeSlab(feasst.args(
#  {"dimension": "2", "bound0": "-1", "bound1": "1"})))))

#mc.add(feasst.Potential(feasst.MakeModelHardShape(feasst.MakeCylinder(
#  feasst.args({"radius": "4"}),
#  feasst.Position(feasst.args({"x": "0", "y": "0", "z": "0"})),
#  feasst.Position(feasst.args({"x": "0", "y": "0", "z": "1"}))
#))))

#mc.add(feasst.Potential(feasst.MakeModelHardShape(feasst.MakeSphere(
#  feasst.args({"radius": "4"}),
#  feasst.Position(feasst.args({"x": "0", "y": "0", "z": "0"}))
#))))

mc.add(feasst.Potential(feasst.MakeModelHardShape(feasst.MakeShapeUnion(
  feasst.MakeSphere(
    feasst.args({"radius": "2.5"}),
    feasst.Position(feasst.args({"x": "0", "y": "0", "z": "0"}))
  ),
  # feasst.MakeSlab(feasst.args({"dimension": "2", "bound0": "-1", "bound1": "1"}))
  feasst.MakeCylinder(
    feasst.args({"radius": "1.5"}),
    feasst.Position(feasst.args({"x": "0", "y": "0", "z": "0"})),
    feasst.Position(feasst.args({"x": "0", "y": "0", "z": "1"}))
  )
))))

mc.add(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": "1.5", "chemical_potential": "1."})))
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "max_move": "2."})))
mc.add(feasst.MakeLog(feasst.args({"steps_per" : "10000"})))
mc.add(feasst.MakeMovie(feasst.args(
  {"steps_per" : "10000", "file_name" : "movie.xyz"})))
mc.add(feasst.MakeEnergyCheck(feasst.args(
  {"steps_per" : "10000", "tolerance" : "1e-10"})))
mc.add(feasst.MakeTuner(feasst.args({"steps_per" : "10000"})))
mc.seek_num_particles(50)
mc.attempt(int(1e6))
