import copy
import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=24.)
parser.add_argument("--dccb_begin", type=int, help="number of molecules before DCCB", default=200)
parser.add_argument("--molecule", type=str, help="type of particle", default=fst.install_dir() + "/forcefield/spce.fstprt")
parser.add_argument("--temperature", type=float, help="temperature in Kelvin", default=525)
parser.add_argument("--max_molecules", type=int, help="maximum number of molecules", default=265)
parser.add_argument("--cubic_box_length", type=float, help="length of cubic periodic boundary conditions", default=20)
parser.add_argument("--beta_mu", type=float, help="beta times chemical potential", default=-8.14)
parser.add_argument("--steps_per", type=int, help="number of MC trials per analysis", default=int(1e5))
args = parser.parse_args()
print("args:", args)

def mc(thread, mn, mx):
    mc = fst.MakeMonteCarlo()
    #mc.set(fst.MakeRandomMT19937(fst.args({"seed": "123"})))
    mc.add(fst.MakeConfiguration(fst.args({"cubic_box_length": str(args.cubic_box_length),
        "physical_constants": "CODATA2010",
        "particle_type0": args.molecule})))
    mc.add(fst.MakePotential(fst.MakeEwald(fst.args({"alpha": str(5.6/args.cubic_box_length),
        "kmax_squared": "38"}))))
    mc.add(fst.MakePotential(fst.MakeModelTwoBodyFactory(fst.MakeLennardJones(), fst.MakeChargeScreened()),
                                      fst.args({"table_size": "1e6"})))
    mc.add(fst.MakePotential(fst.MakeChargeScreenedIntra(), fst.MakeVisitModelBond()))
    mc.add(fst.MakePotential(fst.MakeChargeSelf()))
    mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
    if mx > args.dccb_begin:
        mc.run(fst.MakeAddReference(fst.args({"potential_index": "1", "cutoff": "3.16555789", "use_cell": "true"})))
    beta = 1./fst.kelvin2kJpermol(args.temperature, mc.configuration())
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(beta),
        "chemical_potential": str(args.beta_mu/beta)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        fst.MakeTransitionMatrix(fst.args({"min_sweeps": "10"}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "1.", "tunable_param": "1.",})))
    if mx > args.dccb_begin:
        regrow1 = [{"angle": "true", "mobile_site": "1", "anchor_site": "0", "anchor_site2": "2"}]
        regrow2 = [{"angle": "true", "mobile_site": "2", "anchor_site": "0", "anchor_site2": "1"}]
        regrow12 = [{"bond": "true", "mobile_site": "1", "anchor_site": "0"}] + copy.deepcopy(regrow2)
        regrow21 = [{"bond": "true", "mobile_site": "2", "anchor_site": "0"}] + copy.deepcopy(regrow1)
        grow012 = [{"transfer": "true", "site": "0", "weight": "4"}] + copy.deepcopy(regrow12)
        grow021 = [{"transfer": "true", "site": "0", "weight": "4"}] + copy.deepcopy(regrow21)
        for grow in [regrow1, regrow2]:
            grow[0]["weight"] = "0.3"
        for grow in [regrow12, regrow21]:
            grow[0]["weight"] = "0.2"
        for grow in [grow012, grow021, regrow12, regrow21, regrow1, regrow2]:
            grow[0]["particle_type"] = "0"
            mc.add(fst.MakeTrialGrow(fst.ArgsVector(grow),
                                     fst.args({"reference_index": "0", "num_steps": "4"})))
    else:
        mc.add(fst.MakeTrialRotate(fst.args({"weight": "1.", "tunable_param": "1."})))
        mc.add(fst.MakeTrialTransfer(fst.args({"particle_type": "0", "weight": "4"})))
    mc.add(fst.MakeCheckEnergyAndTune(fst.args({"steps_per": str(args.steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(args.steps_per), "file_name": "clones" + str(thread)})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(args.steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({
        "steps_per": str(args.steps_per),
        "file_name": "clones" + str(thread) + "_crit.txt"})))
    mc.add(fst.MakeEnergy(fst.args({
        "steps_per_write": str(args.steps_per),
        "file_name": "en" + str(thread),
        "multistate": "true"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
    "alpha": "1.75",
    "num": str(args.num_procs),
    "maximum": str(args.max_molecules),
    "extra_overlap": "2"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs)
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt"}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
