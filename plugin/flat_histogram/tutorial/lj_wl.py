import sys
import argparse
import feasst
sys.path.insert(0, '../../core/tutorial') # lj directory
import lj

parser = argparse.ArgumentParser()
parser.add_argument("--box_length", help="periodic cubic box length", default=8, type=float)
parser.add_argument("--temperature", help="temperature", default=1.5, type=float)
parser.add_argument("--chemical_potential", help="chemical potential", default=-2.352321, type=float)
parser.add_argument("--num_steps", help="number of Monte Carlo steps", default=1000000, type=int)
parser.add_argument("--steps_per", help="number of steps per writing log", default=10000, type=int)
args = parser.parse_args()
print("#", args)

def criteria_wl():
  criteria = feasst.MakeCriteriaFlatHistogram(feasst.args(
    {"beta": str(1./args.temperature),
     "chemical_potential": str(args.chemical_potential)}))
  criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
    {"width": "1", "max": "370"}))))
  # criteria.set(feasst.MakeBiasWangLandau(feasst.args({"min_flatness": "20"})))
  criteria.set(feasst.MakeBiasTransitionMatrix(feasst.args(
    {"min_sweeps": "20", "num_steps_to_update": "10000"})))
  return criteria

feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.set(lj.system(box_length=args.box_length))
mc.set(criteria_wl())
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "max_move": "2."})))
mc.add(feasst.MakeTrialTransfer(feasst.args(
  {"weight": "1."})))
lj.add_analysis(mc, args.steps_per)
mc.add(feasst.MakeCriteriaWriter(feasst.args(
  {"steps_per": str(args.steps_per), "file_name": "crit.txt"})))
mc.add(feasst.MakeLog(feasst.args(
  {"steps_per" : str(args.steps_per), "file_name": "log.txt"})))
# mc.seek_num_particles(370)
mc.run_until_complete()
#mc.attempt(args.num_steps)
print(mc.criteria().write())
