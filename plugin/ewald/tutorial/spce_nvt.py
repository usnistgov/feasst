import sys
sys.path.insert(0, '../../monte_carlo/tutorial') # lj directory
import lj
import feasst
import spce

#parser = argparse.ArgumentParser()
#parser.add_argument("--box_length", "-t", help="box_length", default=10., type=float)
#args = parser.parse_args()
#print("#", args)

feasst.seed_random_by_date()
mc = feasst.MonteCarlo()
mc.set(
  spce.system(
    feasst.Configuration(feasst.args(
      {"cubic_box_length": "24.8586887",
       "particle_type": "../../../forcefield/data.spce"})),
    alphaL=5.6,
    kmax_squared=38
  )
)
print("here")
mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
  {"beta": str(1./(298*feasst.ideal_gas_constant/1e3)), # mol/kJ
   "chemical_potential": "0.1"})))
mc.add(feasst.MakeTrialTranslate(feasst.args(
  {"weight": "1.", "tunable_param": "2."})))
mc.add(feasst.MakeTrialRotate(feasst.args(
  {"weight": "1.", "tunable_param": "2."})))
steps_per = int(1e4)
lj.add_analysis(mc, 1)
#lj.add_analysis(mc, steps_per)
#assert(False) # implement random rotation of inserted particles.
mc.seek_num_particles(512)
mc.attempt(int(1e5))
