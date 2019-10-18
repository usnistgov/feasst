import feasst

monte_carlo = feasst.MonteCarlo()
monte_carlo.add(feasst.Configuration(feasst.args(
    {"cubic_box_length": "8",
     "particle_type": feasst.install_dir() + "/forcefield/data.lj"})))
monte_carlo.add(feasst.Potential(feasst.MakeModelLJ()))
monte_carlo.add(feasst.Potential(feasst.MakeLongRangeCorrections()))
monte_carlo.add(feasst.MakeCriteriaMetropolis(feasst.args(
    {"beta": "1.5", "chemical_potential": "1."})))
monte_carlo.add(feasst.MakeTrialTranslate(feasst.args(
    {"tunable_param": "2.", "tunable_target_acceptance": "0.2"})))
monte_carlo.add(feasst.MakeTuner(feasst.args({"steps_per" : "10000"})))
monte_carlo.seek_num_particles(50)
steps_per = int(1e4)
monte_carlo.add(feasst.MakeLog(feasst.args({"steps_per" : str(steps_per)})))
monte_carlo.add(feasst.MakeMovie(feasst.args(
    {"steps_per" : str(steps_per), "file_name" : "movie.xyz"})))
monte_carlo.add(feasst.MakeCheckEnergy(feasst.args(
    {"steps_per" : str(steps_per), "tolerance" : "1e-8"})))
monte_carlo.set(feasst.MakeRandomMT19937(feasst.args({"seed" : "date"})))
monte_carlo.attempt(int(1e6))
