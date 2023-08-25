import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=str, help="random number generator seed",
    default="time")
parser.add_argument("--length", type=float, help="cubic box length", default=8)
parser.add_argument("--num", type=int, help="number of particles", default=50)
parser.add_argument("--data", type=str, help="LMP particle data file",
    default=fst.install_dir() + "/particle/lj.fstprt")
parser.add_argument("--beta", type=float, help="inverse temperature",
    default=1.2)
parser.add_argument("--trials", type=str, help="number of Monte Carlo trials",
    default="1e6")
args = parser.parse_args()
print("args:", args)

# high temperature gcmc to generate initial configuration
mc = fst.MonteCarlo()
mc.set(fst.MakeRandomMT19937(fst.args({"seed" : args.seed})))
mc.add(fst.MakeConfiguration(fst.args({
  "cubic_side_length": str(args.length), "particle_type": args.data})))
mc.add(fst.MakePotential(fst.MakeLennardJones()))
mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
mc.set(fst.MakeThermoParams(fst.args({"beta": "0.1", "chemical_potential": "10"})))
mc.set(fst.MakeMetropolis())
mc.add(fst.MakeTrialTranslate(fst.args(
    {"tunable_param": "2.", "tunable_target_acceptance": "0.2"})))
mc.add(fst.MakeTrialAdd(fst.args({"particle_type": "0"})))
mc.run(fst.MakeRun(fst.args({"until_num_particles": "50"})))

# nvt equilibration
mc.run(fst.MakeRemoveTrial(fst.args({"name": "TrialAdd"})))
mc.set(fst.MakeThermoParams(fst.args({"beta": str(args.beta)})))
mc.add(fst.MakeCheckEnergy(fst.args({"trials_per" : "1e5", "tolerance" : "1e-8"})))
mc.add(fst.MakeTune())
mc.run(fst.MakeRun(fst.args({"num_trials": "1e5"})))

# nvt production
mc.add(fst.MakeLogAndMovie(fst.args({"trials_per" : "1e5", "file_name" : "lj"})))
mc.run(fst.MakeRun(fst.args({"num_trials": args.trials})))
