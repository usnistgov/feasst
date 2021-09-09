import pandas as pd
import matplotlib.pyplot as plt
import feasst as fst
import pyfeasst

# read checkpoints
num_procs = 32
clones = fst.MakeClones('checkpoint', num_procs, 0, ".fst.bak")
volume = clones.clone(0).configuration().domain().volume()
#beta = clones.clone(0).thermo_params().beta()
gce = pyfeasst.find_equilibrium(fst.GrandCanonicalEnsemble(clones), beta_mu_guess=-7.2)


# if checkpoints aren't working:
#volume = 28**3
#beta = 0.306817
#lnpi = pd.read_csv('ln_prob.txt', header=None)
#print(lnpi)
#gce = fst.GrandCanonicalEnsemble(
#    fst.Histogram(fst.args({"width": "1", "max": str(len(lnpi)-1), "min": "0"})),
#    fst.LnProbability(fst.DoubleVector(lnpi[0])),
#    -7)
#gce.reweight(-0.2)
#plt.plot(gce.ln_prob().values())
#plt.show()
#gce = pyfeasst.find_equilibrium(gce, beta_mu_guess=-7.2)

# analyze and plot results
C_MW = 12.0107
H_MW = 1.00784
mw= 4*C_MW + 10*H_MW # C5H12 g/mol
rho_conv = mw/1e6*1e30/fst.CODATA2018().avogadro_constant()
print('<N_vapor>', gce.average_macrostate(0))
print('<rho_vapor>', gce.average_macrostate(0)/volume*rho_conv)
print('<N_liquid>', gce.average_macrostate(1))
print('<rho_vapor>', gce.average_macrostate(1)/volume*rho_conv)
plt.plot(gce.ln_prob().values())
plt.xlabel('N', fontsize=16)
plt.ylabel(r'$\ln\Pi$', fontsize=16)
plt.show()

