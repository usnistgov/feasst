import sys
import argparse
import feasst
sys.path.insert(0, '../../monte_carlo/tutorial') # lj directory
import lj
#import spce
import pyfeasst

parser = argparse.ArgumentParser()
parser.add_argument("--box_length", help="periodic cubic box length", default=8, type=float)
parser.add_argument("--temperature", help="temperature", default=1.5, type=float)
parser.add_argument("--chemical_potential", help="chemical potential", default=-2.352321, type=float)
parser.add_argument("--num_steps", help="number of Monte Carlo steps", default=1000000, type=int)
parser.add_argument("--steps_per", help="number of steps per writing log", default=10000, type=int)
args = parser.parse_args()
print("#", args)

def criteria_flathist(macro_min, macro_max):
  criteria = feasst.MakeCriteriaFlatHistogram(feasst.args(
    {"beta": str(1./args.temperature),
     "chemical_potential": str(args.chemical_potential)}))
  criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
    {"width": "1", "min": str(macro_min), "max": str(macro_max)}))))
  # criteria.set(feasst.MakeBiasWangLandau(feasst.args({"min_flatness": "20"})))
  criteria.set(feasst.MakeBiasTransitionMatrix(feasst.args(
    {"min_sweeps": "20", "num_steps_to_update": "10000"})))
  return criteria

def mc(proc, macro_min=0, macro_max=370, sysflag=1):
  assert(sysflag == 0 or sysflag == 1)
  print("hi")
  mc = feasst.MonteCarlo()
  if sysflag == 0:
    mc.set(lj.system(box_length=args.box_length))
  elif sysflag == 1:
    mc.set(lj.system(box_length=args.box_length))

  # add the minimum number of particles to the system
  mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
    {"beta": "0.001", "chemical_potential": "1"})))
  mc.add(feasst.MakeTrialTranslate(feasst.args(
    {"weight": "1.", "tunable_param": "2."})))
  mc.add(feasst.MakeTrialAdd(feasst.args({"weight": "1."})))
  mc.add(feasst.MakeTrialRemove(feasst.args({"weight": "1."})))
  mc.seek_num_particles(macro_min)

  # replace the acceptance criteria with flat histogram
  mc.set(criteria_flathist(macro_min, macro_max))

  lj.add_analysis(mc, args.steps_per, proc=proc, log="log"+str(proc)+".txt")
  mc.add(feasst.MakeCriteriaWriter(feasst.args(
    {"steps_per": str(args.steps_per), "file_name": "crit"+str(proc)+".txt"})))
  mc.add(feasst.MakeLog(feasst.args(
    {"steps_per" : str(args.steps_per), "file_name": "log"+str(proc)+".txt"})))
  mc.add(feasst.MakeEnergy(feasst.args(
    {"file_name": "energy"+str(proc)+".txt",
     "steps_per_update": "1",
     "steps_per_write": str(args.steps_per),
     "multistate": "true"})))
  mc.add(feasst.MakeCheckpoint(feasst.args(
    {"file_name": "checkpoint"+str(proc)+".txt", "num_hours": "0.01"})))
  mc.run_until_complete()
  #mc.attempt(args.num_steps)
  print(mc.criteria().write())

feasst.seed_random_by_date()

#windows=[[0,2],[2,15],[15,25]]
windows=feasst.window(0, 370, 8, 2)
print(windows)
#windows=pyfeasst.vector_vector_to_list(feasst.window(0, 25, 3, 2))

serial = False
if serial:
  mc(0, 0, 100)
  #for window in windows:
    #mc(proc=0, macro_min=window[0], macro_max=window[1])
  #mc2=mc(macro_min=2, macro_max=15)

else:
  # parallel
  import multiprocessing as mp
  processes = [mp.Process(target=mc, args=(x, windows[x][0], windows[x][1])) for x in range(len(windows))]
  #processes = [mp.Process(target=mc, args=(window[0], window[1])) for window in windows]

  for p in processes:
    p.start()

  for p in processes:
    p.join()

