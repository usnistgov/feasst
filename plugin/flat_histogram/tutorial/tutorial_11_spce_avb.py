import copy
import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=265)
parser.add_argument("--min_particles", type=int, help="minimum number of particles", default=0)
parser.add_argument("--dccb_begin", type=int, help="number of molecules before DCCB", default=200)
parser.add_argument("--avb_end", type=int, help="do not use AVB if the max number of molecules is greater", default=200)
args = parser.parse_args()
print("args:", args)

# determine if avb or dccb is to be used based on maximum number of particles
def avb_or_dccb(mx):
    avb = False
    dccb = False
    if mx > args.dccb_begin:
        dccb = True
    if mx < args.avb_end:
        avb = True
    return avb, dccb

def initialize_neighbor_list(mc):
    neigh_crit = fst.MakeNeighborCriteria(fst.args({"maximum_distance": "10", "minimum_distance": "2.5", "site_type0": "0", "site_type1": "0", "potential_index": "1"}))
    mc.add(neigh_crit)
    lj_and_coul = fst.MakeModelTwoBodyFactory()
    lj_and_coul.add(fst.MakeLennardJones())
    lj_and_coul.add(fst.MakeChargeScreened(fst.args({"table_size": "0"})))
    mc.set(1, fst.MakePotential(lj_and_coul,
        fst.MakeVisitModel(fst.MakeVisitModelInner(fst.MakeEnergyMapNeighborCriteria(neigh_crit))),
        fst.args({"table_size": "1e6"})))

def mc(thread, mn, mx):
    steps_per=int(1e5)
    avb, dccb = avb_or_dccb(mx)
    mc = fst.MakeMonteCarlo()
    spce_args = {"physical_constants": "CODATA2010", "cubic_box_length": "20",
        "alphaL": "5.6", "kmax_squared": "38"}
    if dccb or avb:
        spce_args["dual_cut"] = "3.16555789"
    mc.set(fst.spce(fst.args(spce_args)))
    if avb:
        initialize_neighbor_list(mc)
    beta = 1./fst.kelvin2kJpermol(525, mc.configuration())
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(beta),
        "chemical_potential": str(-8.14/beta)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        fst.MakeTransitionMatrix(fst.args({"min_sweeps": "10"}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "1.", "tunable_param": "1.",})))
    mc.add(fst.MakeTrialRotate(fst.args({"weight": "1.", "tunable_param": "1."})))
    mc.add(fst.MakeTrialTransfer(fst.args({"particle_type": "0", "weight": "4"})))
    regrow1 = [{"angle": "true", "mobile_site": "1", "anchor_site": "0", "anchor_site2": "2"}]
    regrow2 = [{"angle": "true", "mobile_site": "2", "anchor_site": "0", "anchor_site2": "1"}]
    regrow12 = [{"bond": "true", "mobile_site": "1", "anchor_site": "0"}] + copy.deepcopy(regrow2)
    regrow21 = [{"bond": "true", "mobile_site": "2", "anchor_site": "0"}] + copy.deepcopy(regrow1)
    if dccb:
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
    if avb:
        for avbtype in ["transfer_avb", "regrow_avb2", "regrow_avb4"]:
          avb = [{avbtype: "true", "site": "0", "neighbor_index": "0", "target_particle_type": "0", "target_site": "0"}]
          avb_012 = copy.deepcopy(avb) + copy.deepcopy(regrow12)
          avb_021 = copy.deepcopy(avb) + copy.deepcopy(regrow21)
          for grow in [avb_012, avb_021]:
              grow[0]["weight"] = "0.5"
              grow[0]["particle_type"] = "0"
              mc.add(fst.MakeTrialGrow(fst.ArgsVector(grow),
                                       fst.args({"num_steps": "4", "reference_index": "0"})))
    mc.add(fst.MakeCheckEnergyAndTune(fst.args({"steps_per": str(steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeCheckRigidBonds(fst.args({"steps_per": str(steps_per)})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(steps_per), "file_name": "clones" + str(thread)})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({
        "steps_per": str(steps_per),
        "file_name": "clones" + str(thread) + "_crit.txt"})))
    mc.add(fst.MakeEnergy(fst.args({
        "steps_per_write": str(steps_per),
        "file_name": "en" + str(thread),
        "multistate": "true"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
    "alpha": "1.75",
    "num": str(args.num_procs),
    "minimum": str(args.min_particles),
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
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt"}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
