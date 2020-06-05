import feasst as fst

monte_carlo = fst.MonteCarlo()
monte_carlo.set(fst.MakeRandomMT19937(fst.args({"seed" : "time"})))
monte_carlo.add(fst.Configuration(fst.MakeDomain(fst.args({"cubic_box_length": "8"})),
     fst.args({"particle_type": fst.install_dir() + "/forcefield/data.lj"})))
monte_carlo.add(fst.Potential(fst.MakeLennardJones()))
monte_carlo.add(fst.Potential(fst.MakeLongRangeCorrections()))
monte_carlo.add(fst.MakeMetropolis(fst.args({"beta": "1.5"})))
monte_carlo.add(fst.MakeTrialTranslate(fst.args(
    {"tunable_param": "2.", "tunable_target_acceptance": "0.2"})))
steps_per = int(1e3)
monte_carlo.add(fst.MakeTuner(fst.args({"steps_per" : str(steps_per)})))
fst.SeekNumParticles(50)\
    .with_metropolis(fst.args({"beta": "0.1", "chemical_potential": "10"}))\
    .with_trial_add().run(monte_carlo)
monte_carlo.add(fst.MakeLog(fst.args({"steps_per" : str(steps_per)})))
monte_carlo.add(fst.MakeMovie(fst.args(
    {"steps_per" : str(steps_per), "file_name" : "movie.xyz"})))
monte_carlo.add(fst.MakeCheckEnergy(fst.args(
    {"steps_per" : str(steps_per), "tolerance" : "1e-8"})))
monte_carlo.attempt(int(1e5))
