import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--seed", type=str, help="random number generator seed",
    default="time")
parser.add_argument("--length", type=float, help="cubic box length", default=8)
parser.add_argument("--num", type=int, help="number of particles", default=50)
parser.add_argument("--data", type=str, help="LMP forcefield data file",
    default=fst.install_dir() + "/forcefield/data.lj")
parser.add_argument("--beta", type=float, help="inverse temperature",
    default=1.2)
parser.add_argument("--trials", type=int, help="number of Monte Carlo trials",
    default=int(1e6))
args = parser.parse_args()
print("args:", args)

mc = fst.MonteCarlo()
if args.task > 0:
    mc = fst.MakeMonteCarlo("checkpoint.fst")
    mc.attempt(args.trials - mc.trials().num_attempts())
    quit()
mc.set(fst.MakeRandomMT19937(fst.args({"seed" : args.seed})))
mc.add(fst.Configuration(
    fst.MakeDomain(fst.args({"cubic_box_length": str(args.length)})),
    fst.args({"particle_type": args.data})))
mc.add(fst.MakePotential(fst.MakeLennardJones()))
mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
mc.set(fst.MakeThermoParams(fst.args({"beta": str(args.beta)})))
mc.set(fst.MakeMetropolis())
mc.add(fst.MakeTrialTranslate(fst.args(
    {"tunable_param": "2.", "tunable_target_acceptance": "0.2"})))
mc.add(fst.MakeCheckEnergyAndTune(fst.args(
    {"steps_per" : str(int(1e5)), "tolerance" : "1e-8"})))
mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst",
                                    "num_hours": "0.01",
                                    "num_hours_terminate": "0.03"})))
fst.SeekNumParticles(args.num)\
    .with_thermo_params(fst.args({"beta": "0.1", "chemical_potential": "10"}))\
    .with_metropolis()\
    .with_trial_add()\
    .run(mc)
mc.add(fst.MakeLogAndMovie(fst.args(
    {"steps_per" : str(int(1e5)), "file_name" : "lj"})))
mc.attempt(args.trials)
