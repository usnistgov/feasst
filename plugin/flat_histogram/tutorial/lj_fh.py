import sys
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/')
import lj_system
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/tutorial/')
import analyze
sys.path.insert(0, feasst.install_dir() + '/plugin/flat_histogram/tutorial/')
import fh

def monte_carlo(proc=0,                          # processor number
                criteria=fh.criteria_flathist(), # flat histogram criteria
                steps_per=1e6,                   # steps per analysis
                forcefield='data.lj',            # particle type
                run=True                         # run the simulation
                ):
    """Create, run and return a flat histogram grand canonical Monte Carlo simulation"""
    monte_carlo0 = feasst.MonteCarlo()
    monte_carlo0.set(lj_system.system(lj_system.configuration(box_length=8, forcefield=forcefield)))

    # add the minimum number of particles
    monte_carlo0.set(feasst.MakeMetropolis(feasst.args(
        {"beta": "0.001", "chemical_potential": "1"})))
    monte_carlo0.add(feasst.MakeTrialTranslate(feasst.args(
        {"weight": "0.75", "tunable_param": "2."})))
    feasst.add_trial_transfer(monte_carlo0, feasst.args(
        {"weight": "0.125", "particle_type": "0"}))
    monte_carlo0.seek_num_particles(criteria.macrostate().soft_min() + 1)

    # replace the acceptance criteria with flat histogram
    monte_carlo0.set(criteria)

    analyze.add(monte_carlo0, steps_per, proc=proc, log="log"+str(proc)+".txt")

    # periodically write the status of the flat histogram criteria
    monte_carlo0.add(feasst.MakeCriteriaWriter(feasst.args(
        {"steps_per": str(steps_per), "file_name": "crit"+str(proc)+".txt"})))

    # periodically write a checkpoint file
    monte_carlo0.add(feasst.MakeCheckpoint(feasst.args(
        {"file_name": "checkpoint"+str(proc)+".txt", "num_hours": "0.01"})))

    # periodically write the energy of the macrostates
    monte_carlo0.add(feasst.MakeEnergy(feasst.args(
        {"file_name": "energy"+str(proc)+".txt",
         "steps_per_update": "1",
         "steps_per_write": str(steps_per),
         "multistate": "true"})))

    monte_carlo0.run_until_complete()
    #print(monte_carlo0.criteria().write())
    return monte_carlo0
