******************************
Python Interface
******************************

The pybind11 interface is in recent active development.

Installation proceeds as follows:

.. code-block:: bash

    sudo [brew/apt/yum/dnf] install python3-dev pybind11
    # for Rocky, sudo dnf install python3.12-devel.x86_64
    # for Rocky, sudo dnf --enablerepo=devel install python3.12-pybind11
    # or similar
    # replace 3.12 throughout with the relevant version

    # optionally make a virtual environment to juggle multiple versions
    sudo dnf install python3-venv
    mkdir ~/.pyenv
    pushd ~/.pyenv
    python3.12 -m venv feasst

    cd /path/to/feasst/build
    cmake -DUSE_PYBIND11=ON ..
    make -j4
    export LD_LIBRARY_PATH="/path/to/feasst/build/:$LD_LIBRARY_PATH"
    # On Mac, use DYLD_FALLBACK_LIBRARY_PATH
    pip install /path/to/feasst/
    cp _core* ~/.pyenv/feasst/lib/python3.12/site-packages/feasst/
    # The last line didn't use to be necessary.
    # Look for fix in the following links:
    # https://github.com/pybind/scikit_build_example/tree/master
    # https://learn.scientific-python.org/development/
    python ../python/tutorial/test.py

Python usage is as follows:

.. code-block:: py

    import feasst
    mc = feasst.MonteCarlo()
    text_input = """RandomMT19937 seed 123
    Configuration cubic_side_length 8 particle_type0 /feasst/particle/lj.fstprt"""
    for line in text_input.split('\n'):
        feasst.parse(mc, line)
