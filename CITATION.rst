Citation
###############

Hatch HW, Siderius DW, Shen VK (2024) Monte Carlo molecular simulations with FEASST version 0.25.1 *J. Chem. Phys.* 161:092501. https://doi.org/10.1063/5.0224283

Please also site the original papers of the specific methods used.
To aid in this process, a Python tool is in development to recommend citations based on the methods used in the FEASST input text file.
This tool is provided for convenience but does not currently recommend all appropriate citations, and thus the output should not be relied on as complete.
The following command, after replacing $HOME with the appropriate directory of your feasst install, will output recommended citations in bibtex format.

.. code-block:: bash

    python $HOME/feasst/dev/tools/recommend_citations.py --input_file $HOME/feasst/tutorial/example.txt --source_dir $HOME/feasst/
