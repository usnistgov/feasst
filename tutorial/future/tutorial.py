import argparse
import feasst as fst

print(fst.version())
parser = argparse.ArgumentParser()
parser.add_argument("--seed", type=str, help="random number generator seed",
    default="time")
parser.add_argument("--length", type=float, help="cubic box length", default=8)
parser.add_argument("--num", type=int, help="number of particles", default=50)
parser.add_argument("--data", type=str, help="LMP forcefield data file",
    default=fst.install_dir() + "/forcefield/lj.fstprt")
parser.add_argument("--beta", type=float, help="inverse temperature",
    default=1.2)
parser.add_argument("--steps_per", type=str, help="number of trials per analysis or check",
    default="1e5")
parser.add_argument("--equilibration", type=str, help="number of equilibration trials",
    default="1e6")
parser.add_argument("--production", type=str, help="number of production trials",
    default="1e6")
args = parser.parse_args()
print("args:", args)

mc = fst.MakeMonteCarlo(fst.arglist([
    # high temperature gcmc to generate initial configuration
    ["RandomMT19937", {"seed": args.seed}],
    ["Configuration", {"cubic_box_length": str(args.length),
                       "particle_type": args.data}],
    ["Potential", {"Model": "LennardJones"}],
    ["Potential", {"VisitModel": "LongRangeCorrections"}],
    ["ThermoParams", {"beta": "0.1", "chemical_potential": "10"}],
    ["Metropolis", {}],
    ["TrialTranslate", {"tunable_param": "2.", "tunable_target_acceptance": "0.2"}],
    ["TrialAdd", {"particle_type": "0"}],
    ["Run", {"until_num_particles": str(args.num)}],

    # nvt equilibration
    ["RemoveTrial", {"name": "TrialAdd"}],
    ["ThermoParams", {"beta": str(args.beta)}],
    ["Tune", {"steps_per": args.steps_per}],
    ["CheckEnergy", {"steps_per": args.steps_per, "tolerance": "1e-8"}],
    ["Run", {"num_attempts": str(args.equilibration)}],

    # nvt production
    ["RemoveModify", {"name": "Tune"}],
    ["Log", {"steps_per": args.steps_per, "file_name": "lj.txt"}],
    ["Movie", {"steps_per": args.steps_per, "file_name": "lj.xyz"}],
    ["Energy", {"steps_per_write": args.steps_per, "file_name": "en.txt"}],
    ["Run", {"num_attempts": str(args.production)}],
]))
