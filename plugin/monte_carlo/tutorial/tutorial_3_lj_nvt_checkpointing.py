'''
This example, coupled with the slurm script of the same name,
shows an example of using a checkpoint file to restart the simulation
every hour.

Compare the resulting average energy to the NIST SRSW.
https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
'''

import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--seed", type=str, help="random number generator seed",
    default="time")
parser.add_argument("--density", type=float, help="number density, num/volume", default=0.776)
parser.add_argument("--num", type=int, help="number of particles", default=500)
parser.add_argument("--data", type=str, help="LMP forcefield data file",
    default=fst.install_dir() + "/forcefield/lj.fstprt")
parser.add_argument("--beta", type=float, help="inverse temperature",
    default=1/0.9)
parser.add_argument("--equilibration_trials", type=int, help="number of trials for equilibration",
    default=int(5e7))
parser.add_argument("--trials", type=int, help="total number of trials (equilibration and production)",
    default=int(3e8))
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
args = parser.parse_args()
print("args:", args)

mc = fst.MonteCarlo()
if args.task > 0:
    mc = fst.MakeMonteCarlo("checkpoint.fst")
    mc.attempt(args.trials - mc.trials().num_attempts())
    quit()
mc.set(fst.MakeRandomMT19937(fst.args({"seed" : args.seed})))
mc.add(fst.Configuration(
    fst.MakeDomain(fst.args({"cubic_box_length": str((args.num/args.density)**(1./3.))})),
    fst.args({"particle_type": args.data})))
mc.add(fst.MakePotential(fst.MakeLennardJones()))
mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
mc.set(fst.MakeThermoParams(fst.args({"beta": str(args.beta)})))
mc.set(fst.MakeMetropolis())
mc.add(fst.MakeTrialTranslate(fst.args(
    {"tunable_param": "0.2", "tunable_target_acceptance": "0.2"})))
mc.add(fst.MakeTrialAdd(fst.args({"particle_type": "0"})))
mc.run(fst.MakeRun(fst.args({"until_num_particles": str(args.num)})))
mc.run(fst.MakeRemoveTrial(fst.args({"name": "TrialAdd"})))
trials_per = str(int(1e5))
mc.add(fst.MakeCheckEnergy(fst.args({"trials_per" : trials_per, "tolerance" : "1e-8"})))
mc.add(fst.MakeTune())
mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst",
                                    "num_hours": str(0.95*args.num_hours),
                                    "num_hours_terminate": str(0.95*args.num_hours)})))
mc.add(fst.MakeLogAndMovie(fst.args(
    {"trials_per" : trials_per, "file_name" : "lj"})))
mc.add(fst.MakeIncrementPhase(fst.args({"num_trials": str(args.equilibration_trials)})))
mc.add(fst.MakeEnergy(fst.args(
    {"trials_per_write": trials_per, "file_name": "en_lj.txt", "start_after_phase": "0"})))
mc.attempt(args.trials)
