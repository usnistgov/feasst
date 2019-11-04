import feasst

def criteria_flathist(temperature=1.5,
                      chemical_potential=-2.352321,
                      macro_min=0,    # minimum macrostate
                      macro_max=370,  # maximum macrostate
                      tmmc=True,      # use Transition-Matrix (TM) if true, else Wang-Landau (WL)
                      iterations=20,  # number of sweeps (TM) or flatness (WL)
                      steps_per=1e6): # if TM, steps per update
    """Return a flat histogram acceptance criteria with number of particles as the macrostate"""
    criteria = feasst.MakeFlatHistogram(feasst.args(
        {"beta": str(1./temperature),
         "chemical_potential": str(chemical_potential)}))
    criteria.set(feasst.MakeMacrostateNumParticles(feasst.Histogram(feasst.args(
        {"width": "1", "min": str(macro_min), "max": str(macro_max)}))))
    if tmmc:
        criteria.set(feasst.MakeTransitionMatrix(feasst.args(
            {"min_sweeps": str(iterations), "num_steps_to_update": str(steps_per)})))
    else:
        criteria.set(feasst.MakeWangLandau(feasst.args({"min_flatness": str(iterations)})))
    return criteria
