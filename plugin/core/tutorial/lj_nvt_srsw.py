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
  {"weight": "1.", "max_move": "2."})))
lj.add_analysis(mc, args.steps_per)
mc.seek_num_particles(args.num_particles)

# equilibrate
mc.attempt(args.num_steps_equilibrate)

# production
energy = feasst.Accumulator()
energy.set_block(args.steps_per)
for trial in range(args.num_steps):
  mc.attempt(1)
  energy.accumulate(mc.criteria().running_energy())
  if trial % args.steps_per == 0:
    print("energy", energy.average(), "+/-", energy.block_stdev())

print("energy", energy.average(), "+/-", energy.block_stdev())
print("box length", mc.system().configuration().domain().side_length().str())
