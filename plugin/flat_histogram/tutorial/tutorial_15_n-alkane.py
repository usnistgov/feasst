import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--task", type=int, help="SLURM job array index", default=0)
parser.add_argument("--num_procs", type=int, help="number of processors", default=12)
parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=5*24)
parser.add_argument("--cutoff", type=int, help="intermolecular potential cutoff distance", default=14)
parser.add_argument("--dccb_cutoff", type=int, help="cutoff for dual cut configurational bias", default=4)
parser.add_argument("--dccb_begin", type=int, help="number of molecules before using DCCB", default=0)
parser.add_argument("--lx", type=float, help="box length in x", default=28.0)
parser.add_argument("--ly", type=float, help="box length in y", default=28.0)
parser.add_argument("--lz", type=float, help="box length in z", default=28.0)
parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=136)
parser.add_argument("--temperature", type=float, help="temperature in Kelvin", default=392)
parser.add_argument("--particle", "-p", type=str, help="data file of alkane", required=True)
parser.add_argument("--collect_flatness", type=int, help="number of WL flatness to begin collection", default=18)
parser.add_argument("--min_flatness", type=int, help="number of WL flatness to switch to TM", default=22)
parser.add_argument("--beta_mu", type=float, help="baseline chemical potential of each species", default=-7.)
parser.add_argument("--steps_per", type=int, help="steps per analysis", default=int(1e6))
args = parser.parse_args()
print("args:", args)

# return the argument list for growing a complete n-alkane, and in reverse.
def gen_grow(num_sites, reptate=False):
    grow=[]
    grow_inv=[]
    for site in range(num_sites):
        site_inv = num_sites - site - 1
        if reptate:
            if site == num_sites - 1:
                grow += [{"dihedral": "true", "mobile_site": str(site), "anchor_site": str(site-1), "anchor_site2": str(site-2), "anchor_site3": str(site-3)}]
                grow_inv += [{"dihedral": "true", "mobile_site": "0", "anchor_site": "1", "anchor_site2": "2", "anchor_site3": "3"}]
            else:
                grow +=     [{"reptate": "true", "num_steps": "1", "mobile_site": str(site), "anchor_site": str(site+1)}]
                grow_inv += [{"reptate": "true", "num_steps": "1", "mobile_site": str(site_inv), "anchor_site": str(site_inv - 1)}]
        elif site == 0:
            grow += [{"transfer": "true", "site": "0"}]
            grow_inv += [{"transfer": "true", "site": str(site_inv)}]
        elif site == 1:
            grow += [{"bond": "true", "mobile_site": "1", "anchor_site": "0"}]
            grow_inv += [{"bond": "true", "mobile_site": str(site_inv), "anchor_site": str(site_inv + 1)}]
        elif site == 2:
            grow += [{"angle": "true", "mobile_site": "2", "anchor_site": "1", "anchor_site2": "0"}]
            grow_inv += [{"angle": "true", "mobile_site": str(site_inv), "anchor_site": str(site_inv + 1), "anchor_site2": str(site_inv + 2)}]
        else:
            grow += [{"dihedral": "true", "mobile_site": str(site), "anchor_site": str(site-1), "anchor_site2": str(site-2), "anchor_site3": str(site-3)}]
            grow_inv += [{"dihedral": "true", "mobile_site": str(site_inv), "anchor_site": str(site_inv + 1), "anchor_site2": str(site_inv + 2), "anchor_site3": str(site_inv + 3)}]
    return grow, grow_inv

# return a new grow with particle type and weight added to the first arguments
def add_particle_type_weight(site, num_sites, grow):
    import copy
    g= copy.deepcopy(grow[site:num_sites])
    g[0]["particle_type"] = "0"
    weight = 1./float(num_sites - site)
    if site == 0: weight = 4
    g[0]["weight"] = str(weight)
    print(g)
    return g

def mc(thread, mn, mx):
    mc = fst.MakeMonteCarlo()
    mc.add(fst.MakeConfiguration(fst.args({"side_length0": str(args.lx), "side_length1": str(args.ly), "side_length2": str(args.lz),
        "particle_type0": args.particle})))
    for site_type in range(mc.configuration().num_site_types()):
            mc.get_system().get_configuration().set_model_param("cutoff", site_type, args.cutoff)
    mc.add(fst.MakePotential(fst.MakeLennardJones()))
    mc.add(fst.MakePotential(fst.MakeLennardJones(),
                             fst.MakeVisitModelIntra(fst.args({"cutoff": "4"}))))
    mc.add(fst.MakePotential(fst.MakeLongRangeCorrections()))
    if mx > args.dccb_begin:
        reference = fst.Potential(fst.MakeLennardJones(), fst.MakeVisitModelCell(fst.args({"min_length": str(args.dccb_cutoff)})))
        reference.set_model_params(mc.configuration())
        for site_type in range(mc.configuration().num_site_types()):
            reference.set_model_param("cutoff", site_type, args.dccb_cutoff)
        mc.add_to_reference(reference)
        #mc.add_to_reference(fst.MakePotential(fst.MakeLennardJones(),
        #                    fst.MakeVisitModelIntra(fst.args({"cutoff": "4"}))))
        stage_args = {"reference_index": "0", "num_steps": "4"}
    else:
        mc.add_to_reference(fst.MakePotential(fst.DontVisitModel()))
        stage_args = {"reference_index": "0", "num_steps": "1"}
    beta = 1./fst.kelvin2kJpermol(args.temperature);
    mc.set(fst.MakeThermoParams(fst.args({"beta": str(beta), "chemical_potential0": str(args.beta_mu/beta)})))
    mc.set(fst.MakeFlatHistogram(
        fst.MakeMacrostateNumParticles(
            fst.Histogram(fst.args({"width": "1", "max": str(mx), "min": str(mn)}))),
        # fst.MakeTransitionMatrix(fst.args({"min_sweeps": str(args.sweeps)})),
        fst.MakeWLTM(fst.args({
            "collect_flatness": str(args.collect_flatness),
            "min_flatness": str(args.min_flatness),
            "min_sweeps": "1000"}))))
    mc.add(fst.MakeTrialTranslate(fst.args({"weight": "0.5"})))
    mc.add(fst.MakeTrialRotate(fst.args({"weight": "0.5"})))

    # configurational bias with TrialGrow: full and partial regrows from 0-site, and reverse
    num_sites = mc.configuration().particle_type(0).num_sites()
    for site in range(num_sites):
        for g in gen_grow(num_sites, reptate=False):
            mc.add(fst.MakeTrialGrow(fst.ArgsVector(add_particle_type_weight(site, num_sites, g)),
                                     fst.args(stage_args)))

    # reptation
    if num_sites > 3:
        for g in gen_grow(num_sites, reptate=True):
            mc.add(fst.MakeTrialGrow(fst.ArgsVector(add_particle_type_weight(0, num_sites, g)),
                                     fst.args(stage_args)))

    # these moves may not take bond energies into consideration properly. Instead, should implement a dihedral rotation.
    #mc.add(fst.MakeTrialCrankshaft(fst.args(dict({"weight": "0.25", "tunable_param": "25.", "max_length": "5."}, **stage_args))))
    #mc.add(fst.MakeTrialPivot(fst.args(dict({"weight": "0.25", "tunable_param": "25.", "max_length": "5."}, **stage_args))))

    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": str(args.steps_per), "tolerance": "0.0001"})))
    mc.add(fst.MakeTune(fst.args({"steps_per": str(args.steps_per), "stop_after_phase": "0"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": str(args.steps_per),
                                         "file_name": "clones" + str(thread),
                                         "file_name_append_phase": "True"})))
    mc.add(fst.MakeEnergy(fst.args({
        "file_name": "en" + str(thread) + '.txt',
        "file_name_append_phase": "True",
        "start_after_phase": "0",
        "steps_per_write": str(args.steps_per),
        "steps_per_update": "1",
        "multistate": "True"})))
    mc.add(fst.MakeCriteriaUpdater(fst.args({"steps_per": str(args.steps_per)})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"steps_per": str(args.steps_per),
                                            "file_name": "clones" + str(thread) + "_crit.txt",
                                            "file_name_append_phase": "True"})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(thread) + ".fst",
                                        "num_hours_terminate": str(0.9*args.num_procs*args.num_hours)})))
    return mc

windows=fst.WindowExponential(fst.args({
    "alpha": "1.75",
    "num": str(args.num_procs),
    "maximum": str(args.max_particles),
    "extra_overlap": "0"})).boundaries()
print(windows)

if args.task == 0:
    clones = fst.MakeClones()
    for proc, win in enumerate(windows):
        clones.add(mc(proc, win[0], win[1]))
    clones.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint.fst"})))
else:
    clones = fst.MakeClones("checkpoint", args.num_procs)
#clones.initialize_and_run_until_complete()
clones.initialize_and_run_until_complete(fst.args({"ln_prob_file": "ln_prob.txt"}))
print(clones.ln_prob().values())
open('clones.fst', 'w').write(clones.serialize())
