import sys
import math
import unittest
import argparse
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/py') # lj directory
import lj
import pyfeasst

def criteria_flathist(macro_min=0, macro_max=370, chemical_potential=-2.352321, tmmc=True):
    criteria = feasst.MakeCriteriaFlatHistogram(feasst.args(
      {"beta": str(1./1.5),
       "chemical_potential": str(chemical_potential)}))
    criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
      {"width": "1", "min": str(macro_min), "max": str(macro_max)}))))
    if tmmc:
        criteria.set(feasst.MakeBiasTransitionMatrix(feasst.args(
          {"min_sweeps": "20", "num_steps_to_update": str(int(1e6))})))
    else:
        criteria.set(feasst.MakeBiasWangLandau(feasst.args({"min_flatness": "20"})))
    return criteria

def mc(proc=0,
       criteria=criteria_flathist(),
       steps_per=1e6,
       forcefield='data.lj'):
    mc = feasst.MonteCarlo()
    mc.set(lj.system(box_length=8, forcefield=forcefield))

    # add the minimum number of particles to the system
    mc.set(feasst.MakeCriteriaMetropolis(feasst.args(
      {"beta": "0.001", "chemical_potential": "1"})))
    mc.add(feasst.MakeTrialTranslate(feasst.args(
      {"weight": "0.75", "tunable_param": "2."})))
    feasst.add_trial_transfer(mc, feasst.args(
      {"weight": "0.125", "particle_type": "0"}))
    mc.seek_num_particles(criteria.macrostate().soft_min() + 1)

    # replace the acceptance criteria with flat histogram
    mc.set(criteria)

    lj.add_analysis(mc, steps_per, proc=proc, log="log"+str(proc)+".txt")
    mc.add(feasst.MakeCriteriaWriter(feasst.args(
      {"steps_per": str(steps_per), "file_name": "crit"+str(proc)+".txt"})))
    mc.add(feasst.MakeEnergy(feasst.args(
      {"file_name": "energy"+str(proc)+".txt",
       "steps_per_update": "1",
       "steps_per_write": str(steps_per),
       "multistate": "true"})))
    mc.add(feasst.MakeCheckpoint(feasst.args(
      {"file_name": "checkpoint"+str(proc)+".txt", "num_hours": "0.01"})))
    mc.run_until_complete()
    #print(mc.criteria().write())
