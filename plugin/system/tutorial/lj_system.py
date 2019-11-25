"""
This demonstrates the use of modules to re-use various initializations of FEASST
objects in a more pythonic style than the provided C++ interface.
"""

import feasst

def configuration(box_length=8, forcefield='data.lj', num=0):
    """Return an LJ configuration with cubic box and (optionally) cell list.

    box_length -- the length of the cubic peroidic boundary conditions
    forcefield -- the file describing the particle
    num -- the number of particles of the first type to add.
    """
    config = feasst.Configuration(feasst.args(
        {"cubic_box_length": str(box_length),
         "particle_type": feasst.install_dir() + '/forcefield/' + forcefield,
         "init_cells": "3"})) # optionally attempt to create a cell list
    for _ in range(num):
        config.add_particle_of_type(0)
    return config

def system(config):
    """Return an LJ system with long-range corrections and (optionally)
       an optimized potential.

    config -- the configuration
    """
    sys = feasst.System()
    sys.add(config)
    sys.add(feasst.Potential(feasst.MakeLennardJones()))
    sys.add(feasst.Potential(feasst.MakeLongRangeCorrections()))
    # If the requested cell list spacing is not large enough for the box, the cell list will
    # be disabled even if a cell list was requested previously.
    if sys.configuration().domain().is_cell_enabled():
        # While adding a potential is by default considered an "unoptimized"
        # potential, a system may also contain an optimized potential.
        # Optimized potentials are used for most energy calculations, but
        # CheckEnergy asserts that it is equivalent to the unoptimized
        # potential.
        sys.add_to_optimized(feasst.Potential(feasst.MakeLennardJones(), feasst.MakeVisitModelCell()))
        sys.add_to_optimized(feasst.Potential(feasst.MakeLongRangeCorrections()))
        # While the cell list is optimized, the LRCs are not. Regardless, the LRCs are added so
        # that the energy is equivalent to the unoptimized potentials.
    return sys
