import sys
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/')
import lj_system
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/tutorial/')
import analyze
sys.path.insert(0, feasst.install_dir() + '/plugin/flat_histogram/tutorial/')
import fh

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--nopara", help="disable prefetch", dest='nopara', action='store_true')
parser.set_defaults(nopara=False)
parser.add_argument("--nofh", help="disable fh", dest='nofh', action='store_true')
parser.set_defaults(nofh=False)
parser.add_argument("--nobalance", help="disable load_balance", dest='nobalance', action='store_true')
parser.set_defaults(nobalance=False)
parser.add_argument("--cutoff", type=float, default=4.)
parser.add_argument("--target_prob", type=float, default=0.2)
parser.add_argument("--rel_disp_prob", type=float, default=5)
parser.add_argument("--max_move", type=float, default=0.185)
parser.add_argument("--density", type=float, default=0.85)
parser.add_argument("--temperature", type=float, default=0.88)
parser.add_argument("--chemical_potential", type=float, default=-2.35)
parser.add_argument("--window_half_width", type=int, default=10)
parser.add_argument("--iterations", type=int, default=20)
args = parser.parse_args()
print(args)

print(feasst.install_dir())
assert(feasst.install_dir() == '/home/hwh/gcfetch/feasst')

def lj_system(box_length):
    system = feasst.System()
    config = feasst.Configuration(
        feasst.Domain(feasst.args({"cubic_box_length": str(box_length)})),
        feasst.args({"particle_type": feasst.install_dir() + "/forcefield/lj.fstprt"}))
    config.set_model_param("cutoff", 0, args.cutoff)
    system.add(config)
    system.add(feasst.Potential(feasst.MakeLennardJones()))
    return system

def mc(
       steps_per=int(1e6),
       ):
    box_length = 2.*args.cutoff
    file_app = "_a" + str(args.rel_disp_prob) + "_rc" + str(args.cutoff)

    monte_carlo = feasst.Prefetch(feasst.args({"steps_per_check": str(int(1e7))}))
    monte_carlo.activate_prefetch(False)
    # monte_carlo.set(feasst.MakeRandomMT19937(feasst.args({"seed": "1578687129"})))
    monte_carlo.set(lj_system(box_length=box_length))
    monte_carlo.set(feasst.MakeMetropolis(feasst.args({
        "beta": str(1./args.temperature),
        "chemical_potential": str(args.chemical_potential),
    })))
    monte_carlo.add(feasst.MakeTrialTranslate(feasst.args({
        "weight": str(args.rel_disp_prob),
        "tunable_param": str(args.max_move),
        "tunable_target_acceptance": str(args.target_prob),
        "tunable_percent_change": "0.1",
    })))
    num_particles = int(args.density*box_length**3)
    nmin = num_particles - args.window_half_width
    nmax = num_particles + args.window_half_width
    if not args.nofh:
        feasst.SeekNumParticles(nmin).with_trial_add().run(monte_carlo)
    monte_carlo.add(feasst.MakeTrialTransfer(feasst.args({
        "weight": "1",
        "particle_type": "0"})))
    if not args.nofh:
        monte_carlo.set(fh.criteria_flathist(
            temperature=args.temperature,
            chemical_potential=args.chemical_potential,
            macro_max=nmax,
            macro_min=nmin,
            iterations=args.iterations,
            ))
        monte_carlo.add(feasst.MakeCriteriaUpdater(feasst.args({"steps_per": str(steps_per)})))
        monte_carlo.add(feasst.MakeCriteriaWriter(feasst.args(
            {"steps_per": str(steps_per), "file_name": "crit"+file_app+".txt"})))
    else:
        monte_carlo.add(feasst.MakeNumParticles(feasst.args({
            "file_name": "num"+file_app+".txt",
            "steps_per_write": str(steps_per),
        })))
    analyze.add(monte_carlo,
        steps_per,
        proc=file_app,
        log="log"+file_app+".txt",
        )

    if not args.nopara:
        monte_carlo.activate_prefetch(True)

    monte_carlo.add(feasst.MakeCPUTime(feasst.args({
        "steps_per_update": str(steps_per),
        "steps_per_write": str(steps_per),
        "file_name": "cpu" + file_app + ".txt",
    })))

    monte_carlo.add(feasst.MakeEnergy(feasst.args(
        {"file_name": "energy"+file_app+".txt",
         "steps_per_update": "1",
         "steps_per_write": str(steps_per),
         "multistate": "true"})))

    monte_carlo.set(feasst.MakeCheckpoint(feasst.args(
        {"file_name": "checkpoint"+file_app+".txt", "num_hours": "0.1"})))

    # run until complete is not pprefetched correctly
    monte_carlo.run_until_complete()
    #monte_carlo.attempt(int(1e7))
mc()
