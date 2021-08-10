import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=512)
parser.add_argument("--cubic_box_length", type=float, help="side length of cubic periodic boundaries", default=8)
parser.add_argument("--mu", type=float, help="chemical potential", default=-2.352321)
parser.add_argument("--min_sweeps", type=int, help="minimum number of TM sweeps before termination", default=100)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    steps_per=str(int(1e6))
    mc = fst.MakeMonteCarlo()
    mc.add(fst.Configuration(
      fst.MakeDomain(fst.args({"cubic_box_length": str(args.cubic_box_length)})),
      fst.args({"particle_type0": fst.install_dir() + "/forcefield/data.hard_sphere"})))
    mc.add(fst.MakePotential(fst.MakeHardSphere()))
#    mc.add_to_optimized(fst.MakePotential(fst.MakeHardSphere(), fst.MakeVisitModelCell(fst.args({"min_length": "1"}))))
    mc.set(fst.MakeThermoParams(fst.args({"beta": "1",
                                          "chemical_potential": str(args.mu)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.min_sweeps)})),
        fst.MakeWLTM(fst.args({"collect_flatness": "18",
                               "min_flatness": "22",
                               "min_sweeps": str(args.min_sweeps)}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"new_only": "true", "weight": "1.", "tunable_param": "1."})))
    mc.add(fst.MakeTrialTransfer(fst.args({"new_only": "true", "weight": "4", "particle_type": "0"})))
    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": steps_per, "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"steps_per": steps_per, "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": steps_per,
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakePairDistribution(fst.args({"steps_per_update": "1000",
                                              "steps_per_write": steps_per,
                                              "dr": "0.025",
                                              "file_name": "gr" + str(thread) + ".txt",
                                              "file_name_append_phase": "True",
                                              "start_after_phase": "0",
                                              "multistate": "True",
                                              "multistate_aggregate": "False"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": steps_per})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"steps_per": steps_per,
                                            "file_name": "clones" + str(thread) + "_crit.txt",
                                            "file_name_append_phase": "True"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours": str(0.1*args.num_procs*args.num_hours),
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
    "alpha": "2.5",
    "num": str(args.num_procs),
    "maximum": str(args.max_particles)})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs)
#clones.initialize_and_run_until_complete()
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt",
                                                   "omp_batch": str(int(1e6))}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
