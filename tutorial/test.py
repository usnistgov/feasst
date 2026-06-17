import numpy as np
import pandas as pd
import feasst

# Compare with T*=0.9,rho*=0.003 in https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.
num_particles=500
density=0.003
beta=1./0.9
length=np.power(num_particles/density, 1./3.)
prefix='lj'
write=f'trials_per_write=1e4 output_file={prefix}'

mc = feasst.MonteCarlo()
feasst.parse(mc, f"""RandomMT19937 seed=416974832
Configuration cubic_side_length={length} particle_type=lj:/feasst/particle/lj_new.txt
Potential Model=LennardJones VisitModel=VisitModelCell
ThermoParams beta={beta} chemical_potential=1
Metropolis
TrialTranslate tunable_param=2
Checkpoint checkpoint_file={prefix}_checkpoint.fst num_hours=1 num_hours_terminate=117.5667
CheckEnergy trials_per_update=1e4 decimal_places=8
Log {write}_eq.csv
Movie {write}_eq.xyz
Run until_num_particles={num_particles} Trial=TrialAdd particle_type=lj
Metropolis trials_per_cycle=1e4 cycles_to_complete=10
Run until=complete Stepper=Tune
Remove name=Log,Movie""")

print('# x-position of first particle/site.', mc.configuration(0).particle(0).site(0).position(0))
assert mc.configuration(0).particle(0).site(0).position(0) != 0.

feasst.parse(mc, f"""Metropolis trials_per_cycle=1e4 cycles_to_complete=1e2
Log {write}.csv
Movie {write}.xyz
Energy {write}_en.csv
CPUTime {write}_cpu.csv
ProfileCPU {write}_profile.csv
GhostTrialVolume {write}_pressure.csv trials_per_update=1e4
Run until=complete""")

print('# compare average energy with SRSW https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm.')
en = pd.read_csv('lj_en.csv')
print('<U> =', en['average'][0]/num_particles, '+/-', en['block_stdev'][0]/num_particles)
