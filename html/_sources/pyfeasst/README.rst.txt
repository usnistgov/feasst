pyfeasst
===========================================================

Python package to run and analyze FEASST simulations.
Install pyfeasst locally to match your local FEASST version:

.. code-block:: bash

    pip install /path/to/feasst/pyfeasst

Juggle multiple pyfeasst versions with Python virtual environments:

.. code-block:: bash

    sudo [apt,yum,brew,dnf] install python3-dev python3-venv
    mkdir ~/.pyenv
    cd ~/.pyenv
    python3 -m venv feasst
    source ~/.pyenv/feasst/bin/activate # may add this to your .bash_profile

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   docs/source/accumulator
   docs/source/cd
   docs/source/fstio
   docs/source/fstplot
   docs/source/macrostate_distribution
   docs/source/multistate_accumulator
   docs/source/physical_constants

