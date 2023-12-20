#from pyfeasst import physical_constants
r_0_OO = 3.5532 # A
eps_OO = 0.1848 # kcal/mol
A = eps_OO*r_0_OO**12
print('A_OO (kcal A^12/mol)', A)
eps_HH = 0.01 # kcal/mol
j_per_cal = 4.184
print('eps_OO (kJ/mol)', eps_OO*j_per_cal)
print('eps_HH (kJ/mol)', eps_HH*j_per_cal)
k_HOH = 60 # kcal/mol/rad^2
print('k_HOH (kJ/mol/rad^2)', k_HOH*j_per_cal)
k_OH = 250 # kcal/mol/A^2
print('k_OH (kJ/mol/A^2)', k_OH*j_per_cal)
