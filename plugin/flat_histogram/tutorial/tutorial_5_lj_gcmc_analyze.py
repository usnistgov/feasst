import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import feasst as fst
import pyfeasst

num_procs = 12
clones = fst.MakeClones('checkpoint', num_procs)
beta = clones.clone(0).thermo_params().beta()
volume = clones.clone(0).configuration().domain().volume()
plt.plot(clones.ln_prob().values(), label='T*='+str(1./beta))
sat=list()
for temp in np.arange(0.8, 1.201, 0.05):
    extrap = fst.ExtrapolateBetaGCE(clones,
        fst.args({"beta_new": str(1/temp),
                  "beta_original": str(beta)}))
    extrap=pyfeasst.find_equilibrium(extrap, beta_mu_guess=-6)
    plt.plot(extrap.ln_prob().values(), label='T*='+str(round(temp,2)))

    # tabulate saturation properties
    num_vapor = extrap.average_macrostate(0)
    num_liquid = extrap.average_macrostate(1)
    sat.append([temp,
        num_vapor/volume,
        num_liquid/volume,
        extrap.betaPV()/volume/beta,
        extrap.average(extrap.energy(), 0)/num_vapor,
        extrap.average(extrap.energy(), 1)/num_liquid])

print(pd.DataFrame(sat, columns=['T', 'rhov', 'rhol', 'pressure', 'uv', 'ul']))
plt.legend()
plt.xlabel(r'$N$', fontsize=16)
plt.ylabel(r'$\ln\Pi$', fontsize=16)
plt.show()
