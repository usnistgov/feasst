.. Edit /feasst/dev/sphinx/README.rst.in (CMake template) and not /feasst/README.rst

*************************
README
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free and open-source software to conduct molecular- and particle-based simulations with Monte Carlo methods.
New users can start with the `website <https://pages.nist.gov/feasst/>`_ (`DOI <https://doi.org/10.18434/M3S095>`_), `manuscript <https://doi.org/10.1063/5.0224283>`_, `GitHub discussion <https://github.com/usnistgov/feasst/discussions>`_ and a `five minute video <https://www.nist.gov/video/how-use-feasst-0255-monte-carlo-molecular-simulation-software>`_.
Support FEASST with a `GitHub <https://github.com/usnistgov/feasst>`_ star or `manuscript <https://doi.org/10.1063/5.0224283>`_ citation!

Quickstart
===========

.. code-block:: bash

    sudo apt install g++ cmake python3-venv               # C++, CMake and Python3
    python3 -m venv feasst; source ~/feasst/bin/activate  # or use existing env/conda
    CMAKE_BUILD_PARALLEL_LEVEL=8 pip3 install feasst      # install
    feasst-menu                                           # interactive tutorial

For Mac/Linux use apt, yum, dnf, brew(homebrew), etc.
For Windows, use WSL.
For HPC, use "module avail/load."

C++ Executable
================================

Run with BASH/etc:

.. code-block:: bash

    #!/bin/bash
    num_particles=500
    density=0.003
    beta=`python3 -c "print(1./0.9)"`
    length=`python3 -c "print(($num_particles/$density)**(1./3.))"`
    tpc=1e4
    feasst << EOF
    # Comments begin with the # symbol.
    # Compute average energy of LJ at T*=0.9, rho*=0.003
    # See https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
    MonteCarlo
    RandomMT19937
    #seed=1572362164
    Configuration cubic_side_length=$length particle_type=lj:/feasst/particle/lj.txt
    Potential Model=LennardJones VisitModel=VisitModelCell
    Potential VisitModel=LongRangeCorrections
    Checkpoint checkpoint_file=checkpoint.fst
    ThermoParams beta=$beta chemical_potential=-1
    Metropolis
    TrialTranslate weight=1 tunable_param=2
    Tune
    CheckEnergy trials_per_update=$tpc decimal_places=8
    Log trials_per_write=$tpc output_file=lj_eq.csv format=vertical
    Run until_num_particles=$num_particles particle_type=lj Trial=TrialAdd weight=2
    Run num_trials=1e5
    Remove name=Tune,Log
    WriteCheckpoint
    EOF
    
Restart with checkpoint.fst and BASH/etc:

.. code-block:: bash

    #!/bin/bash
    write="trials_per_write=1e4 output_file=lj"
    feasst << EOF
    Restart checkpoint_file=checkpoint.fst
    Log $write.csv
    Movie $write.xyz
    Energy ${write}_en.csv
    Metropolis trials_per_cycle=1e4 cycles_to_complete=1e2
    Run until=complete
    EOF
    
Run with script.txt and BASH/etc:

.. code-block:: bash

    feasst < script.txt

Python Module
===============================

Run in Python:

.. code-block:: python

    import numpy as np
    import pandas as pd
    import feasst
    
    # Compare with T*=0.9,rho*=0.003 in https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.
    num_particles=500
    density=0.003
    beta=1./0.9
    length=np.power(num_particles/density, 1./3.)
    prefix='lj'
    write=f'trials_per_write=1e4 output_file={prefix}'
    
    mc = feasst.MonteCarlo()
    feasst.parse(mc, f"""RandomMT19937 seed=416974832
    Configuration cubic_side_length={length} particle_type=lj:/feasst/particle/lj_new.txt
    Potential Model=LennardJones VisitModel=VisitModelCell
    ThermoParams beta={beta} chemical_potential=1
    Metropolis
    TrialTranslate tunable_param=2
    Checkpoint checkpoint_file={prefix}_checkpoint.fst num_hours=1 num_hours_terminate=117.5667
    CheckEnergy trials_per_update=1e4 decimal_places=8
    Log {write}_eq.csv
    Movie {write}_eq.xyz
    Run until_num_particles={num_particles} Trial=TrialAdd particle_type=lj
    Metropolis trials_per_cycle=1e4 cycles_to_complete=10
    Run until=complete Stepper=Tune
    Remove name=Log,Movie""")
    
    print('# x-position of first particle/site.', mc.configuration(0).particle(0).site(0).position(0))
    assert mc.configuration(0).particle(0).site(0).position(0) != 0.
    
    feasst.parse(mc, f"""Metropolis trials_per_cycle=1e4 cycles_to_complete=1e2
    Log {write}.csv
    Movie {write}.xyz
    Energy {write}_en.csv
    CPUTime {write}_cpu.csv
    ProfileCPU {write}_profile.csv
    GhostTrialVolume {write}_pressure.csv trials_per_update=1e4
    Run until=complete""")
    
    print('# compare average energy with SRSW https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.')
    en = pd.read_csv('lj_en.csv')
    print('<U> =', en['average'][0]/num_particles, '+/-', en['block_stdev'][0]/num_particles)
    
Documentation
=======================

See https://pages.nist.gov/feasst. View documentation for your specific version using feasst-menu, or by `downloading FEASST <https://github.com/usnistgov/feasst/tags>`_ and opening the html/index.html file as shown by the following Bash commands.

.. code-block:: bash

    curl -OL https://github.com/usnistgov/feasst/archive/refs/tags/v0.25.19.tar.gz
    tar -xf v0.25.19.tar.gz
    open feasst-0.25.19/html/index.html
