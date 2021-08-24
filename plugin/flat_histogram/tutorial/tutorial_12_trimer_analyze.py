import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import feasst as fst
import pyfeasst

num_procs = 4
clones = fst.MakeClones('checkpoint', num_procs)
beta = clones.clone(0).thermo_params().beta()
volume = clones.clone(0).configuration().domain().volume()
fh = clones.flat_histogram(0)
#plt.plot(clones.ln_prob().values(), label='T*='+str(1./beta))
#plt.show()
sat=list()
for temp in [0.3, 0.29, 0.28, 0.27, 0.26]:
#for temp in np.arange(0.18, 0.3, 0.025):
    print('temp', temp)
    extrap = fst.ExtrapolateBetaGCE(clones,
        fst.args({"beta_new": str(1/temp),
                  "beta_original": str(beta)}))
    num_smooth = 60
    if temp == 0.3: num_smooth = 10
    extrap=pyfeasst.find_equilibrium(extrap, beta_mu_guess=-1, num_smooth=num_smooth)
    plt.plot(extrap.ln_prob().values(), label='T*='+str(round(temp,2)))

# need to pass num_smooth to average_macrostate
#    # tabulate saturation properties
#    num_vapor = extrap.average_macrostate(0)
#    num_liquid = extrap.average_macrostate(1)
#    sat.append([temp,
#        num_vapor/volume,
#        num_liquid/volume,
#        extrap.betaPV()/volume/beta,
#        extrap.average(extrap.energy(), 0)/num_vapor,
#        extrap.average(extrap.energy(), 1)/num_liquid])

print(pd.DataFrame(sat, columns=['T', 'rhov', 'rhol', 'pressure', 'uv', 'ul']))
plt.legend()
plt.xlabel(r'$N$', fontsize=16)
plt.ylabel(r'$\ln\Pi$', fontsize=16)
plt.show()
