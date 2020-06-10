import feasst as fst

mc = fst.MonteCarlo()
mc.set(fst.MakeRandomMT19937(fst.args({"seed" : "time"})))
mc.add(fst.Configuration(fst.MakeDomain(fst.args({"cubic_box_length": "8"})),
     fst.args({"particle_type": fst.install_dir() + "/forcefield/data.lj"})))
mc.add(fst.Potential(fst.MakeLennardJones()))
mc.add(fst.Potential(fst.MakeLongRangeCorrections()))
mc.add(fst.MakeMetropolis(fst.args({"beta": "1.5"})))
mc.add(fst.MakeTrialTranslate(fst.args(
    {"tunable_param": "2.", "tunable_target_acceptance": "0.2"})))
steps_per = int(1e3)
mc.add(fst.MakeCheckEnergyAndTune(fst.args(
    {"steps_per" : str(steps_per), "tolerance" : "1e-8"})))
fst.SeekNumParticles(50)\
    .with_metropolis(fst.args({"beta": "0.1", "chemical_potential": "10"}))\
    .with_trial_add().run(mc)
mc.add(fst.MakeLogAndMovie(fst.args(
    {"steps_per" : str(steps_per), "file_name" : "lj"})))
mc.attempt(int(1e5))
