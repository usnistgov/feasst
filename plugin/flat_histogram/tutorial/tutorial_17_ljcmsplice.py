import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=370)
parser.add_argument("--temperature", type=float, help="temperature", default=1.5)
parser.add_argument("--mu", type=float, help="chemical potential", default=-2.352321)
parser.add_argument("--min_sweeps", type=int, help="minimum number of TM sweeps before termination", default=100)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx, soft_min, soft_max):
    steps_per=str(int(1e6))
    mc = fst.MakeMonteCarlo()
    mc.add(fst.MakeConfiguration(fst.args({"cubic_box_length": "8",
        "particle_type0": fst.install_dir() + "/forcefield/lj.fstprt"})))
    mc.add(fst.MakePotential(fst.MakeLennardJones()))
    mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
    trial_args = {"particle_type": "0", "site": "0"}
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(1./args.temperature),
                                          "chemical_potential": str(args.mu)})))
    mc.set(fst.MakeMetropolis())
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([dict({"translate": "true", "tunable_param": "1"}, **trial_args)])))
    mc.add(fst.MakeTrialAdd(fst.args({"particle_type": "0"})))
    mc.run(fst.MakeRun(fst.args({"until_num_particles": str(soft_min)})))
    mc.run(fst.MakeRemoveTrial(fst.args({"name": "TrialAdd"})))
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(1./args.temperature),
                                          "chemical_potential": str(args.mu)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)})),
            fst.args({"soft_max": str(soft_max), "soft_min": str(soft_min)})),
        fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.min_sweeps)}))))
#        fst.MakeWLTM(fst.args({"collect_flatness": "18",
#                               "min_flatness": "22",
#                               "min_sweeps": str(args.min_sweeps)}))))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([dict({"transfer": "true", "weight": "4"}, **trial_args)])))
    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": steps_per, "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"steps_per": steps_per, "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": steps_per,
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({"steps_per_write": steps_per,
                                    "file_name": "en" + str(thread) + ".txt",
                                    "file_name_append_phase": "True",
                                    "start_after_phase": "0",
                                    "multistate": "True"})))
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
    "maximum": str(args.max_particles),
    "overlap": "0"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeCollectionMatrixSplice()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, 0, args.max_particles, win[0], win[1]))
    #clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs)
#clones.initialize_and_run_until_complete()
clones.run(0.01)
while not clones.are_all_complete():
    clones.adjust_bounds(15)
    clones.run(0.01)
    open('lnpi.txt', 'w').write(str(clones.ln_prob().values()))
    data=str()
    for i in range(clones.num()):
        data += str(clones.clone(i).criteria().num_iterations()) + " "
        data += str(clones.flat_histogram(i).macrostate().soft_min()) + " "
        data += str(clones.flat_histogram(i).macrostate().soft_max()) + " "
    open('meta.txt', 'a').write(data+"\n")
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
