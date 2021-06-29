import sys
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/')
import lj_system
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/tutorial/')
import analyze

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--nopipe", help="disable pipeline", dest='nopipe', action='store_true')
parser.set_defaults(nopipe=False)
parser.add_argument("--cutoff", type=float)
parser.add_argument("--target_prob", type=float)
parser.add_argument("--max_move", type=float)
args = parser.parse_args()
print(args)

# assert(feasst.install_dir() == '/home/local/NIST/hwh/feasst1.0')
# assert(feasst.install_dir() == '/home/local/NIST/hwh/msd2/feasst')
assert(feasst.install_dir() == '/home/hwh/msd4/feasst')

#accept_to_param={
#  #0.75: 0.04,
#  #0.5: 0.085,
#  #0.45: 0.099,
#  #0.4: 0.111,
#  0.35: 0.125,
#  #0.3: 0.14,
#  0.25: 0.16,
#  #0.225: 0.1725,
#  #0.2: 0.185,
#  #0.175: 0.2,
#  0.15: 0.215,
#  #0.125: 0.2375,
#  #0.1: 0.26,
#  }

#res = dict()

def lj_system(box_length, cutoff):
    system = feasst.System()
    config = feasst.Configuration(
        feasst.MakeDomain(feasst.args({"cubic_box_length": str(box_length)})),
        feasst.args({"particle_type": feasst.install_dir() + '/forcefield/data.lj'}))
    config.set_model_param("cutoff", 0, cutoff)
    system.add(config)
    system.add(feasst.Potential(feasst.MakeLennardJones()))
    if system.configuration().domain().is_cell_enabled():
        #system.add_to_optimized(feasst.Potential(feasst.MakeLennardJones(), feasst.MakeVisitModelCell()))
        system.add_to_reference(feasst.Potential(feasst.MakeHardSphere(), feasst.MakeVisitModelCell()))
    return system

def mc(target_acceptance=0.25,
       tunable_param=0.1,
       density=0.85,
       cutoff=4,
       temperature=0.88,
       steps_per=int(1e6),
       ):
    box_length = 2*cutoff
    file_app = "_a" + str(target_acceptance) + "_r" + str(cutoff)
    num_particles=int(density*box_length**3)
    print('num_particles', num_particles)
    print('target_acceptance', target_acceptance)
    print('tunable_param', tunable_param)

    monte_carlo = feasst.Prefetch(feasst.args({"steps_per_check": "10000000"}))
    monte_carlo.activate_prefetch(False)
    monte_carlo.set(lj_system(box_length=box_length, cutoff=cutoff))
    monte_carlo.set(feasst.MakeMetropolis(feasst.args({
        "beta": str(1./temperature),
        "chemical_potential": "1.",
    })))
    monte_carlo.add(feasst.MakeTrialTranslate(feasst.args({
        "weight": "1.",
        "tunable_param": str(tunable_param),
        "tunable_target_acceptance": str(target_acceptance),
        "tunable_percent_change": "0.01",
        # "num_steps": "4",
        # "reference_index": "0",
    })))
    feasst.SeekNumParticles(num_particles).with_trial_add().run(monte_carlo)
    monte_carlo.add(feasst.MakeLog(feasst.args(
        {"steps_per" : str(steps_per),
         "file_name": "log"+file_app+".txt",
         "clear_file": "true"})))
    monte_carlo.add(feasst.MakeCheckEnergy(feasst.args(
        {"steps_per" : str(steps_per),
         "tolerance" : str(1e-8)})))
    monte_carlo.add(feasst.MakeTune(feasst.args(
        {"steps_per" : str(steps_per)})))

    #equilibrate
    monte_carlo.attempt(int(1e7))

    if not args.nopipe:
        monte_carlo.activate_prefetch(True)
    monte_carlo.add(feasst.MakeMeanSquaredDisplacement(feasst.args({
        "steps_per_update": "10000",
        "updates_per_origin": "1000",
        "file_name": "msd" + file_app + ".txt",
        "steps_per_write": str(int(1e5))
    })))

    monte_carlo.add(feasst.MakeCPUTime(feasst.args({
        "steps_per_update": str(steps_per),
        "steps_per_write": str(steps_per),
        "file_name": "cpu" + file_app + ".txt",
    })))

    #h0 = feasst.cpu_hours()
    monte_carlo.attempt(int(1e8))
    #print("elapsed hours:", feasst.cpu_hours() - h0)
    #res[target_acceptance]=dict()
    #res[target_acceptance]["time"] = feasst.cpu_hours() - h0

#for cutoff in [5]:
#    for acc in accept_to_param:
mc(target_acceptance=args.target_prob,
   tunable_param=args.max_move,
   cutoff=args.cutoff
   )
