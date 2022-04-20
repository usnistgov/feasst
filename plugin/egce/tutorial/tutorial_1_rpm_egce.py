import math
import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--dccb_begin", type=int, help="begin DCCB at this many particles", default=400)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=650)
parser.add_argument("--temperature", type=float, help="temperature", default=0.047899460618081)
parser.add_argument("--beta_mu", type=float, help="chemical potential", default=-13.94)
parser.add_argument("--min_sweeps", type=int, help="minimum number of TM sweeps before termination", default=100)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    trials_per=int(1e4)
    mc = fst.MakeMonteCarlo()
    mc.add(fst.MakeConfiguration(fst.args({"cubic_box_length": "12",
        "particle_type0": fst.install_dir() + "/plugin/charge/forcefield/rpm_plus.fstprt",
        "particle_type1": fst.install_dir() + "/plugin/charge/forcefield/rpm_minus.fstprt",
        "cutoff": "4.891304347826090",
        "charge0": str( 1./math.sqrt(fst.CODATA2018().charge_conversion())),
        "charge1": str(-1./math.sqrt(fst.CODATA2018().charge_conversion()))})))
    mc.add(fst.MakePotential(fst.MakeEwald(fst.args({"alpha": str(6.87098396396261/12),
        "kmax_squared": "38"}))))
    mc.add(fst.MakePotential(fst.MakeModelTwoBodyFactory(fst.MakeHardSphere(), fst.MakeChargeScreened()),
                                      fst.args({"table_size": "1e6"})))
    mc.add(fst.MakePotential(fst.MakeChargeSelf()))
    num_steps = "1"
    ref = "-1"
    if mx >= args.dccb_begin:
        mc.run(fst.MakeAddReference(fst.args({"potential_index": "1", "cutoff": "1", "use_cell": "true"})))
        ref = "0"
        num_steps = "4"
    beta = 1./args.temperature
    mc.set(fst.MakeThermoParams(fst.args({
        "beta": str(beta),
        "chemical_potential0": str(args.beta_mu/beta),
        "chemical_potential1": str(args.beta_mu/beta)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.min_sweeps)})),
        fst.MakeWLTM(fst.args({"collect_flatness": "18",
                               "min_flatness": "22",
                               "min_sweeps": str(args.min_sweeps)})),
        fst.MakeAEqualB(fst.args({"extra_A": "1"}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "1.", "tunable_param": "1."})))
    for particle_type in range(mc.configuration().num_particle_types()):
        mc.add(fst.MakeTrialTransfer(fst.args({"weight": "4", "particle_type": str(particle_type),
            "reference_index": ref, "num_steps": num_steps})))
    mc.add(fst.MakeCheckEnergy(fst.args({"trials_per": str(trials_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"trials_per": str(trials_per), "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"trials_per": str(trials_per),
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({"trials_per_write": str(trials_per),
                                    "file_name": "en" + str(thread) + ".txt.",
                                    "file_name_append_phase": "True",
                                    "start_after_phase": "0",
                                    "multistate": "True"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"trials_per": str(trials_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"trials_per": str(trials_per),
                                            "file_name": "clones" + str(thread) + "_crit.txt",
                                            "file_name_append_phase": "True"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours": str(0.1*args.num_procs*args.num_hours),
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
  "alpha": "2",
  "num": str(args.num_procs),
  "maximum": str(args.max_particles),
  "extra_overlap": "2"})).boundaries()
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
                                                   "omp_batch": str(int(1e4))}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
