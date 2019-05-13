import argparse
import feasst
import lj

parser = argparse.ArgumentParser()
parser.add_argument("--box_length", help="periodic cubic box length", default=8, type=float)
parser.add_argument("--temperature", help="temperature", default=1.5, type=float)
parser.add_argument("--chemical_potential", help="chemical potential", default=-8.352321, type=float)
parser.add_argument("--num_steps", help="number of Monte Carlo steps", default=1000000, type=int)
parser.add_argument("--steps_per", help="number of steps per writing log", default=100000, type=int)
args = parser.parse_args()
print("#", args)

feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.set(lj.system(box_length=args.box_length))
mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": str(1./args.temperature),
   "chemical_potential": str(args.chemical_potential)})))
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "max_move": "2."})))
mc.add(feasst.MakeTrialTransfer(feasst.args(
  {"weight": "1."})))
lj.add_analysis(mc, args.steps_per)
mc.attempt(args.num_steps)
