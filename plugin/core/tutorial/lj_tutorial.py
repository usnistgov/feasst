
import argparse
import feasst

parser = argparse.ArgumentParser()
parser.add_argument("--box_length", "-l", help="box length", default=0, type=float)
parser.add_argument("--num_particles", help="number of particles", default=0, type=int)
parser.add_argument("--density", "-r", help="density. Overrides box length or num particles", default=0., type=float)
parser.add_argument("--temperature", "-t", help="temperature", default=1.5, type=float)
parser.add_argument("--num_steps_equilibrate", help="number of equilibration Monte Carlo steps", default=0, type=int)
parser.add_argument("--num_steps", help="number of Monte Carlo steps", default=0, type=int)
parser.add_argument("--steps_per", help="number of steps per writing log", default=0, type=int)
args = parser.parse_args()
print("#", args)

if args.density != 0.:
  assert(args.box_length == 0.)
  args.box_length = (args.num_particles/args.density)**(1./3.)

def lj():
  potential = feasst.Potential()
  potential.set_model(feasst.ModelLJShrPtr())
  return potential

def lj_cell():
  potential = lj()
  potential.set_visit_model(feasst.VisitModelCell())
  return potential

def lrc():
  potential = feasst.Potential()
  potential.set_visit_model(feasst.LongRangeCorrectionsShrPtr())
  return potential

def system():
  system = feasst.System()
  system.add(feasst.Configuration(feasst.args(
    {"cubic_box_length": str(args.box_length),
     "particle_type": "../../../forcefield/data.lj",
     "init_cells": "3"})))
  system.add(lj())
  system.add_to_optimized(lj_cell())
  system.add(lrc())
  system.add_to_optimized(lrc())
  return system

def monte_carlo():
  mc = feasst.MonteCarlo()
  mc.set(system())
  mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
    {"beta": str(1./args.temperature), "add_activity": "1."})))
  mc.add(feasst.MakeTrialTranslate(feasst.args(
    {"weight": "1.", "max_move": "2."})))
  mc.add(feasst.MakeLog(feasst.args(
    {"steps_per" : str(args.steps_per)})))
  mc.add(feasst.MakeMovie(feasst.args(
    {"steps_per" : str(args.steps_per),
     "file_name" : "movie.xyz"})));
  mc.add(feasst.MakeEnergyCheck(feasst.args(
    {"steps_per" : str(args.steps_per),
     "tolerance" : "1e-10"})));
  mc.add(feasst.MakeTuner(feasst.args(
    {"steps_per" : str(args.steps_per)})));
  return mc

feasst.seed_random_by_date()
mc = monte_carlo()
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
