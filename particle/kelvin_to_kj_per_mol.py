from pyfeasst import physical_constants
epsilon_in_Kelvin = 28.129 # emp2 co2 C-C
epsilon_in_Kelvin = 47.588 # emp2 co2 C-O
epsilon_in_Kelvin = 80.507 # emp2 co2 O-O
print('epsilon(K)', epsilon_in_Kelvin)
print('epsilon(kJ/mol)', epsilon_in_Kelvin*physical_constants.MolarGasConstant().value()/1e3)
