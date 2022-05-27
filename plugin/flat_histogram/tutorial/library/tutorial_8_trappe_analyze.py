import feasst as fst
import pyfeasst

# read checkpoints
num_procs = 12
clones = fst.MakeClones('checkpoint', num_procs)
volume = clones.clone(0).configuration().domain().volume()
beta = clones.clone(0).thermo_params().beta()
gce = pyfeasst.find_equilibrium(fst.GrandCanonicalEnsemble(clones))

# compute saturation pressure
R=clones.clone(0).configuration().physical_constants().ideal_gas_constant()
na=clones.clone(0).configuration().physical_constants().avogadro_constant()
press_conv=R/1e3*1e30/na
print('saturation pressure (kPa)', gce.betaPV()/volume/beta*press_conv)

# compute saturation compositions
num0 = fst.DoubleVector()
clones.stitch(num0, "NumParticles", fst.AccumulatorAverage())
num_vapor = gce.average(num0, 0)
num_liquid = gce.average(num0, 1)

print('vapor y_C2H4', 1 - num_vapor/gce.average_macrostate(0))
print('liquid x_C2H4', 1 - num_liquid/gce.average_macrostate(1))

# obtain extensive moments
extmom_index = fst.SeekAnalyze().index("ExtensiveMoments", clones.clone(0))[0]

# Return ExtensiveMoments, dervied class of Analyze, by serialization, which is relatively slow.
# see steppers/include/ExtensiveMoments.h for moments API
def extensive_moment(window, state):
    return fst.ExtensiveMoments(clones.clone(window).analyze(extmom_index).analyze(state))

extmom = extensive_moment(2, 30)
for p in range(3):
    for m in range(2):
        for k in range(2):
            for j in range(1):
                for i in range(1):
                    print(p, m, k, j, i,
                          extmom.moments(p, m, k, j, i).num_values(),
                          extmom.moments(p, m, k, j, i).sum_dble(),
                          extmom.moments(p, m, k, j, i).sum_of_squared_dble())
