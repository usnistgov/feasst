import sys
import feasst
sys.path.insert(0, feasst.install_dir() + '/plugin/system/tutorial/')
import lj_system
sys.path.insert(0, feasst.install_dir() + '/plugin/monte_carlo/tutorial/')
import analyze

def criteria_flathist(temperature=1.5,
                      chemical_potential=-2.352321,
                      macro_min=0,    # minimum macrostate
                      macro_max=370,  # maximum macrostate
                      tmmc=True,      # use Transition-Matrix (TM) if true, else Wang-Landau (WL)
                      iterations=20): # number of sweeps (TM) or flatness (WL)
    """Return a flat histogram acceptance criteria with number of particles as the macrostate"""
    criteria = feasst.MakeFlatHistogram(feasst.args(
        {"beta": str(1./temperature),
         "chemical_potential": str(chemical_potential)}))
    criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
        {"width": "1", "min": str(macro_min), "max": str(macro_max)}))))
    if tmmc:
        criteria.set(feasst.MakeTransitionMatrix(feasst.args({"min_sweeps": str(iterations)})))
    else:
        criteria.set(feasst.MakeWangLandau(feasst.args({"min_flatness": str(iterations)})))
    return criteria

def monte_carlo(proc=0,                          # processor number
                criteria=criteria_flathist(), # flat histogram criteria
                steps_per=1e5,                   # steps per analysis
                system=lj_system.system(lj_system.configuration(box_length=8, forcefield="data.lj")),
                run=True                         # run the simulation
                ):
    """Create, run and return a flat histogram grand canonical Monte Carlo simulation"""
    monte_carlo0 = feasst.MonteCarlo()
    monte_carlo0.set(system)

    # add the minimum number of particles
    monte_carlo0.set(feasst.MakeMetropolis(feasst.args(
        {"beta": "0.001", "chemical_potential": "1"})))
    monte_carlo0.add(feasst.MakeTrialTranslate(feasst.args(
        {"weight": "0.375", "tunable_param": "2."})))
    if monte_carlo0.system().configuration().particle_type(0).num_sites() > 1:
        monte_carlo0.add(feasst.MakeTrialRotate(feasst.args(
            {"weight": "0.375", "tunable_param": "2."})))
    feasst.add_trial_transfer(monte_carlo0, feasst.args(
        {"weight": "0.125", "particle_type": "0"}))

    monte_carlo0.seek_num_particles(criteria.macrostate().soft_min() + 1)

    # replace the acceptance criteria with flat histogram
    monte_carlo0.set(criteria)

    analyze.add(monte_carlo0, steps_per, proc=proc, log="log"+str(proc)+".txt")

    # periodically write the status of the flat histogram criteria
    monte_carlo0.add(feasst.MakeCriteriaUpdater(feasst.args({"steps_per": str(steps_per)})))
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

    if run:
        monte_carlo0.run_until_complete()
    #print(monte_carlo0.criteria().write())
    return monte_carlo0
