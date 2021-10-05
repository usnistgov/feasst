'''
The NVT+W method is a simplification of Transition Matrix Monte Carlo where each "window" is only one macrostate.
This means that no insertions or deletions are accepted, although they are still attempted in order to compute the transition probabilities.
While NVT+W is more easy to parallelize, it is not recommended because it has major drawbacks.
NVT+W is less efficient because the attempted insertion and deletion moves are never accepted.
This reduces sampling of the Markov chain.
In addition, NVT+W is not a flat-histogram method.
'''
import argparse
from itertools import repeat
from multiprocessing import Pool
import feasst as fst

def nvtw(num_particles, num_procs, num_equil, num_prod, num_hours, dccb_begin, temperature, mu, steps_per):
    mc = fst.MakeMonteCarlo()
    #mc.set(fst.MakeRandomMT19937(fst.args({"seed": "1633373856"})))
    sys_args = dict()
    num_steps = "1"
    ref = "-1"
    if num_particles >= dccb_begin:
        sys_args["dual_cut"] = str(1)
        ref = "0"
        num_steps = "4"
    #mc.set(fst.lennard_jones(fst.args(sys_args)))
    config = fst.MakeConfiguration(fst.args({"cubic_box_length": "8", "particle_type0": fst.install_dir() + "/forcefield/atom.fstprt"}))
    config.set_model_param("cutoff", 0, 1.5)
    mc.add(config)
    mc.add(fst.MakePotential(fst.MakeSquareWell()))

    mc.set(fst.MakeThermoParams(fst.args({"beta": str(1./temperature),
                                          "chemical_potential": str(mu)})))
    mc.set(fst.MakeMetropolis());
    trial_args = {"particle_type": "0", "site": "0", "reference_index": ref, "num_steps": num_steps}
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([dict({"translate": "true", "tunable_param": "1"}, **trial_args)])))
    mc.add(fst.MakeTrialAdd(fst.args({"particle_type": "0", "weight": "4"})))
    mc.add(fst.MakeTune(fst.args({"steps_per": steps_per})))
    mc.add(fst.MakeCheckEnergy(fst.args({"steps_per": steps_per, "tolerance": "0.0001"})))
    mc.add(fst.MakeLogAndMovie(fst.args({"steps_per": steps_per,
                                         "file_name": "lj" + str(num_particles)})))
    mc.set(fst.MakeCheckpoint(fst.args({"file_name": "checkpoint" + str(num_particles) + ".fst",
                                        "num_hours": str(0.1*num_procs*num_hours),
                                        "num_hours_terminate": str(0.9*num_procs*num_hours)})))
    mc.perform(fst.MakeRun(fst.args({"until_num_particles": str(num_particles)})))
    mc.perform(fst.MakeRemoveTrial(fst.args({"name": "TrialAdd"})))
    mc.attempt(int((num_particles+1)*num_equil))
    mc.perform(fst.MakeRemoveModify(fst.args({"name": "Tune"})))
    mc.add(fst.MakeTrialGrow(fst.ArgsVector([dict({"transfer": "true", "weight": "4"}, **trial_args)])))
    mc.set(fst.MakeFlatHistogram(fst.args({
        "Macrostate": "MacrostateNumParticles", "width": "1", "max": str(num_particles), "min": str(num_particles),
        "Bias": "TransitionMatrix", "min_sweeps": "1"})))
    mc.add(fst.MakeEnergy(fst.args({"steps_per_write": steps_per,
                                    "file_name": "en" + str(num_particles) + ".txt"})))
    mc.add(fst.MakeCriteriaWriter(fst.args({"steps_per": steps_per,
                                            "file_name": "crit" + str(num_particles) + ".txt"})))
    mc.attempt(int((num_particles+1)*num_prod))

def run_parallel(num_procs, num_equil, num_prod, num_hours, dccb_begin, max_particles, temperature, mu, steps_per):
    #nums = [max_particles]
    #nvtw(nums[0], num_procs, num_equil, num_prod, num_hours, dccb_begin, temperature, mu, steps_per)
    nums = reversed(range(0, max_particles + 1))
    with Pool(num_procs) as pool:
        pool.starmap(nvtw, zip(nums, repeat(num_procs), repeat(num_equil), repeat(num_prod), repeat(num_hours), repeat(dccb_begin), repeat(temperature), repeat(mu), repeat(steps_per)))

if __name__ == '__main__':
    print(fst.version())
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_procs", type=int, help="number of processors", default=1)
    parser.add_argument("--num_equil", type=int, help="number of equilibration steps per particle", default=1e3)
    parser.add_argument("--num_prod", type=int, help="number of production steps per particle", default=1e4)
    parser.add_argument("--num_hours", type=float, help="number of hours before restart", default=1.)
    parser.add_argument("--dccb_begin", type=int, help="begin DCCB at this many particles", default=3000)
    parser.add_argument("--max_particles", type=int, help="maximum number of particles", default=370)
    parser.add_argument("--temperature", type=float, help="temperature", default=1.5)
    parser.add_argument("--mu", type=float, help="chemical potential", default=-2.352321)
    parser.add_argument("--steps_per", type=str, help="number of trials per analysis", default="1e4")
    args = parser.parse_args()
    print("args:", args)
    run_parallel(**vars(args))

