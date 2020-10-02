import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--num_steps", type=int, help="number of DCCB steps", default=1)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    steps_per=int(1e4)
    mc = fst.MakeMonteCarlo()
    ref = -1
    lj_args = dict()
    if args.num_steps > 1:
        lj_args["dual_cut"] = str(1)
        ref = 0
    mc.set(fst.lennard_jones(fst.args(lj_args)))
    mc.add(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": "10"})),
        fst.MakeWLTM(fst.args({"collect_flatness": "18",
                               "min_flatness": "22",
                               "min_sweeps": "10"})),
        fst.args({"beta": str(1./1.5), "chemical_potential": "-2.352321"})))
    mc.add(fst.MakeTrialTranslate(fst.args({
        "weight": "1.",
        "tunable_param": "1.",
        "reference_index": str(ref),
        "num_steps": str(args.num_steps)})))
    mc.add(fst.MakeTrialTransfer(fst.args({
        "particle_type": "0",
        "weight": "4",
        "reference_index": str(ref),
        "num_steps": str(args.num_steps)})))
    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": str(steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeTuner(fst.args({"steps_per": str(steps_per), "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(steps_per),
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({"steps_per_write": str(steps_per),
                                    "file_name": "en" + str(thread) + ".txt.",
                                    "file_name_append_phase": "True",
                                    "start_after_phase": "0",
                                    "multistate": "True"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"steps_per": str(steps_per),
                                            "file_name": "clones" + str(thread) + "_crit.txt",
                                            "file_name_append_phase": "True"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours": str(0.1*args.num_procs*args.num_hours),
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
  "alpha": "2.25",
  "num": str(args.num_procs),
  "maximum": "370",
  "extra_overlap": "2"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs);
#clones.initialize_and_run_until_complete()
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt",
                                                   "omp_batch": str(int(1e4))}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
