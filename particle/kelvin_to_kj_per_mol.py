from pyfeasst import physical_constants
#epsilon_in_Kelvin = 98 # CH3 2-methylpropane
#epsilon_in_Kelvin = 10 # CH 2-methylpropane
epsilon_in_Kelvin = 27.0 # TraPPE co2 C-C
#epsilon_in_Kelvin = 28.129 # emp2 co2 C-C
#epsilon_in_Kelvin = 47.588 # emp2 co2 C-O
#epsilon_in_Kelvin = 80.507 # emp2 co2 O-O
print('epsilon(K)', epsilon_in_Kelvin)
print('epsilon(kJ/mol)', epsilon_in_Kelvin*physical_constants.MolarGasConstant().value()/1e3)

k_in_Kelvin_per_rad_sq = 0.5*62500 # 2-methylpropane, note FEASST doesn't have the k/2 prefactor, but a k prefactor
print('k(K/rad_sq)', k_in_Kelvin_per_rad_sq)
print('k(kJ/mol/rad_sq)', k_in_Kelvin_per_rad_sq*physical_constants.MolarGasConstant().value()/1e3)

