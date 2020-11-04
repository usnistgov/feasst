import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--dccb_begin", type=int, help="number of molecules before DCCB", default=300)
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    steps_per=int(1e5)
    mc = fst.MakeMonteCarlo()
    spce_args = {"physical_constants": "CODATA2010", "cubic_box_length": "20",
        "alphaL": "5.6", "kmax_squared": "38"}
    if mx > args.dccb_begin: spce_args["dual_cut"] = "3.16555789"
    mc.set(fst.spce(fst.args(spce_args)))
    beta = 1./fst.kelvin2kJpermol(525, mc.configuration())
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(beta),
        "chemical_potential": str(-8.14/beta)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        fst.MakeTransitionMatrix(fst.args({"min_sweeps": "10"}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "1.", "tunable_param": "1.",})))
    if mx > args.dccb_begin:
        ref = "0"
        num = "4"
        mc.add(fst.MakeTrialGrow(fst.ArgsVector([
            {"transfer": "true",  # "regrow": "true",  # regrow isn't very efficient
             "particle_type": "0", "site": "0", "weight": "4"},
            {"bond": "true", "mobile_site": "1", "anchor_site": "0"},
            {"angle": "true", "mobile_site": "2", "anchor_site": "0", "anchor_site2": "1"}]),
            fst.args({"reference_index": ref, "num_steps": num})))
        mc.add(fst.MakeTrialGrow(fst.ArgsVector([
            {"particle_type": "0", "weight": "0.2",
             "bond": "true", "mobile_site": "1", "anchor_site": "0"},
            {"angle": "true", "mobile_site": "2", "anchor_site": "0", "anchor_site2": "1"}]),
            fst.args({"reference_index": ref, "num_steps": num})))
        mc.add(fst.MakeTrialGrow(fst.ArgsVector([
            {"particle_type": "0", "weight": "0.3", "reference_index": ref, "num_steps": num,
             "angle": "true", "mobile_site": "2", "anchor_site": "0", "anchor_site2": "1",}])))
        mc.add(fst.MakeTrialGrow(fst.ArgsVector([
            {"particle_type": "0", "weight": "0.3", "reference_index": ref, "num_steps": num,
             "angle": "true", "mobile_site": "1", "anchor_site": "0", "anchor_site2": "2"}])))
        mc.add(fst.MakeRecenterParticles(fst.args({"steps_per": str(steps_per)})))
        mc.add(fst.MakeCheckRigidBonds(fst.args({"steps_per": str(steps_per)})))
    else:
        mc.add(fst.MakeTrialRotate(fst.args({"weight": "1.", "tunable_param": "1."})))
        mc.add(fst.MakeTrialTransfer(fst.args({"particle_type": "0", "weight": "4"})))
    mc.add(fst.MakeCheckEnergyAndTune(fst.args({"steps_per": str(steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(steps_per), "file_name": "clones" + str(thread)})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({
        "steps_per": str(steps_per),
        "file_name": "clones" + str(thread) + "_crit.txt"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
  "alpha": "1.75",
  "num": str(args.num_procs),
  "maximum": "265",
  "extra_overlap": "2"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs);
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt"}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
