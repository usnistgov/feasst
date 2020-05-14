#!/usr/bin/env python

import feasst
import pyfeasst
from scipy import stats

from FEASST_Henry_Coefficient_Rigid import *

#Relevant CODATA Constants
amu = 1.6605390400e-27 #kg/amu
kB = 1.380649030e-23 #J/K
Na = 6.022140760E23  #1/mol
h = 6.6260701500E-34  #J s

#Confidence Level
conf_level = 0.95

#Temperature
T = 300.

#Reference State Conditions
p = 1. # bar
V = kB * T / (p * 1.e5)  #molecular volume in m3


test_input = { "adsorptive": "data.C2",
               "adsorbent": "./LTA_replicate",
               "temperature": T,
               "rcut": 15.,
               "ncoeffs": 3,
               "trials": 1.e6,
               "pair_type": "LJCoulEwald",
               "tail_type": "LFS",
               "Ewald":{ "k2max": 27, "alpha": 6.0},
               "scale": {"active": True, "factor": 10000.}
             }


# Run the Henry-law Calculations
#  Serial Version
#nthreads = 1
#Kcoeffs, Kcoeffs_var, beta = HenryCoefficient(debug=False,seed=835279,**test_input)
#  Parallel Version
nthreads = 4
Kcoeffs, Kcoeffs_var, beta = Parallel_HenryCoefficient(nthreads=nthreads,seed=835279,**test_input)

# Compute the cell mass and volume from the XYZ file
# NOTE: This is specific to SiO2 cells.
cell_mass = 0.
with open(test_input["adsorbent"]+'.xyz',mode='r') as f:
    lines = f.readlines()
    for i,line in enumerate(lines):
        if i == 1:
             cell = [float(x) for x in line.split()]
        elif i > 1:
            atom = line.split()[0]
            if atom == 'Si':
                cell_mass += 28.0855  #amu
            elif atom == 'O':
                cell_mass += 15.999   #amu
cell_volume = cell[0]*cell[1]*cell[2] #ang^3
# print('mass:', cell_mass)
# print('vol: ', cell_volume)
rhoS = (cell_mass * amu) / (cell_volume * 1.e-30) #Skeletal Density in kg/m3
# print('density:', rhoS)

# Calculate the Henry Coefficient
Kh = Kcoeffs[0] * beta / rhoS
# Units:
#   beta = mol/kJ
#   rhoS = kg/m3
#   Kcoeffs[0] = dimensionless
#   Kh = (mol/kJ)/(kg/m3) = (mol/kg)(m3/1000J) = (mmol/kg)(1/Pa)
#     that is, millimoles adsorbate per kilogram adsorbent per Pascal
# Uncertainty in Kh
nu_Ki = int(test_input["trials"] * nthreads)  #degrees of freedom  in the MC integrals
a0 = (beta/rhoS)
var_Kh = (a0**2)*Kcoeffs_var[0]
coverage_factor = stats.t.ppf(1-(1.-conf_level)/2., nu_Ki-1)
CI_Kh = np.sqrt(var_Kh/float(nu_Ki-1)) * coverage_factor

# Calculate the infinite-dilution Isosteric Heat
Qst_inf = 1./beta + Kcoeffs[1]/Kcoeffs[0]  #FEASST native units = kJ/mol
#   Uncertainty in Qst_inf
a0 = Kcoeffs[1]/(Kcoeffs[0]**2)
a1 = 1./Kcoeffs[0]
var_Qst_inf = (a0**2)*Kcoeffs_var[0] + (a1**2)*Kcoeffs_var[1]

nu_eff = (var_Qst_inf**2)/(
    (a0**4)*(Kcoeffs_var[0]**2)/float(nu_Ki) + (a1**4)*(Kcoeffs_var[1]**2)/float(nu_Ki)
    )
nu_eff = int(nu_eff)
coverage_factor = stats.t.ppf(1-(1.-conf_level)/2., nu_eff-1)
CI_Qst_inf = np.sqrt(var_Qst_inf/float(nu_eff-1)) * coverage_factor

#NOTE: Uncertainty estimates are based on linear uncertainty propagation. Degrees of
# freedom are estimated from the sample count for direct measurements or using
# the Welch-Satterwaithe equation for derived quantities
# Uncertainty is assumed to derive exclusively from the statistical variation
# in the Monte Carlo averages. I.E., no error arises from the values of beta,
# rhoS, etc.

# Output to screen
print(test_input["adsorptive"])
print(test_input["adsorbent"])
print('   Henry Coefficient            = '+str(Kh)+' mmol/(kg Pa)')
print('     +/-                        = '+str(CI_Kh)+' mmol/(kg Pa)')
print('   Isosteric Heat (inf. dilute) = '+str(Qst_inf)+' kJ/mol')
print('     +/-                        = '+str(CI_Qst_inf)+' kJ/mol')
print(' +/- ranges are the 95 % confidence intervals estimated via linear uncertainty propagation')

# Assert that the resulting values are within expected range
assert(abs(Kh - 0.17017) < 3*CI_Kh)
assert(abs(Qst_inf - 33.7794) < 3*CI_Qst_inf)
