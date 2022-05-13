import math
import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
#parser.add_argument("--dccb_begin", type=int, help="begin DCCB at this many particles", default=300)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=370)
parser.add_argument("--temperature", type=float, help="temperature", default=0.7)
parser.add_argument("--mu", type=float, help="chemical potential", default=-1.5)
parser.add_argument("--min_sweeps", type=int, help="minimum number of TM sweeps before termination", default=100)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    trials_per=str(int(1e5))
    mc = fst.MakeMonteCarlo()
    #mc.set(fst.MakeRandomMT19937(fst.args({"seed": "12345"})))
    #mc.set(fst.MakeRandomMT19937(fst.args({"seed": "1634926410"})))
    chi = 0.7 # patch coverage
    patch_angle = 2*math.asin(math.sqrt(chi/2))*180/math.pi
    print('patch_angle', patch_angle)
    mc = fst.MonteCarlo()
    config = fst.MakeConfiguration(fst.args({"cubic_box_length": "8",
        "particle_type0": fst.install_dir() + "/plugin/patch/forcefield/two_patch_linear.fstprt"}))
    config.add(fst.MakeGroup(fst.args({"site_type0": "0"})))
    mc.add(config)
    mc.add(fst.MakePotential(fst.MakeHardSphere(),
                             fst.MakeVisitModelCell(fst.args({"min_length": "1", "cell_group": "1"})),
                             fst.args({"group_index": "1"})))
    patch = fst.MakeVisitModelInnerPatch(fst.args({"patch_degrees_of_type1": str(patch_angle)}))
    mc.add(fst.MakePotential(
        fst.MakeSquareWell(),
        fst.MakeVisitModel(patch),
        #fst.MakeVisitModelCell(patch, fst.args({"min_length": "1.5", "cell_group": "1"})),
        fst.args({"group_index": "1"})))
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(1./args.temperature),
                                          "chemical_potential": str(args.mu)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.min_sweeps)})),
        fst.MakeWLTM(fst.args({"collect_flatness": "18",
                               "min_flatness": "22",
                               "min_sweeps": str(args.min_sweeps)}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "0.5", "tunable_param": "0.1"})))
    mc.add(fst.MakeTrialRotate(fst.args({"weight": "0.5", "tunable_param": "20."})))
    mc.add(fst.MakeTrialTransfer(fst.args({"weight": "4", "particle_type": "0"})))
#    mc.add(fst.MakeTrialGrow(fst.ArgsVector([{"translate": "true", "site": "0", "tunable_param": "1"}]),
#        {"reference_index": ref, "num_steps": num_steps}))
#    mc.add(fst.MakeTrialGrow(fst.ArgsVector([{"transfer": "true", "site": "0", "weight": "4"}]),
#        {"reference_index": ref, "num_steps": num_steps}))
    mc.add(fst.MakeCheckEnergy(fst.args({"trials_per": trials_per, "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"stop_after_phase": "0"})))
    mc.add(fst.MakeLog(fst.args({"trials_per": trials_per,
                                 "file_name": "clones" + str(thread) + ".txt",
                                 "file_name_append_phase": "True"})))
    mc.add(fst.MakeMoviePatch(fst.args({"trials_per": trials_per,
                                        "file_name": "clones" + str(thread) + ".xyz",
                                        "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({"trials_per_write": trials_per,
                                    "file_name": "en" + str(thread) + ".txt",
                                    "file_name_append_phase": "True",
                                    "start_after_phase": "0",
                                    "multistate": "True"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"trials_per": trials_per})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"trials_per": trials_per,
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
clones.initialize_and_run_until_complete()
#clones.initialize(1)
#print("num 0", clones.clone(0).configuration().num_particles())
#print("num 1", clones.clone(1).configuration().num_particles())
#print("current 0", clones.clone(0).criteria().current_energy())
#print("current 1", clones.clone(1).criteria().current_energy())
#clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt",
#                                                   "omp_batch": str(int(1e6))}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
