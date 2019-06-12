import argparse
import feasst
import lj

parser = argparse.ArgumentParser()
parser.add_argument("--num_particles", help="number of particles", default=500, type=int)
parser.add_argument("--density", "-r", help="density. Overrides box length or num particles", default=0.001, type=float)
parser.add_argument("--temperature", "-t", help="temperature", default=0.9, type=float)
parser.add_argument("--num_steps_equilibrate", help="number of equilibration Monte Carlo steps", default=50000000, type=int)
parser.add_argument("--num_steps", help="number of Monte Carlo steps", default=200000000, type=int)
parser.add_argument("--steps_per", help="number of steps per writing log", default=100000, type=int)
args = parser.parse_args()
print("#", args)

feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.set(lj.system(box_length=(args.num_particles/args.density)**(1./3.)))
mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": str(1./args.temperature),
   "chemical_potential": "1."})))
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "tunable_param": "2."})))
lj.add_analysis(mc, args.steps_per)
mc.seek_num_particles(args.num_particles)

# equilibrate
mc.attempt(args.num_steps_equilibrate)

# compute average energy using a stepper/analysis
energy = feasst.MakeEnergy(feasst.args({
  "steps_per_update": "1",
  "steps_per_write": str(args.steps_per),
#  "file_name": "energy.txt",
}))
mc.add(energy);

# compute average using just this script
energy_alt = feasst.Accumulator()

# production
for trial in range(args.num_steps):
  mc.attempt(1)
  energy_alt.accumulate(mc.criteria().current_energy())

print("energy:", energy.energy().str(), "alternative:", energy_alt.average(), "+/-", energy_alt.block_stdev())
print("box length", mc.system().configuration().domain().side_length().str())
