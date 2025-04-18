# From SI Table 1 of https://doi.org/10.1038/s41467-018-08222-6
from pyfeasst import physical_constants
#epsilon_in_eV = 0.268376 # mW
epsilon_in_eV = 0.297284 # ML-mW
print('epsilon (eV)', epsilon_in_eV)
# (eV/molecule) (e J/eV) (kJ/1e3/J) (na molecules/mol)
print('epsilon (kJ/mol)', epsilon_in_eV*physical_constants.ElementaryCharge().value()*physical_constants.AvogadroConstant().value()/1e3)
