#!/bin/bash
num_particles=500
density=0.003
beta=`python3 -c "print(1./0.9)"`
length=`python3 -c "print(($num_particles/$density)**(1./3.))"`
tpc=1e4
feasst << EOF
# Comments begin with the # symbol.
# Compute average energy of LJ at T*=0.9, rho*=0.003
# See https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
MonteCarlo
RandomMT19937
#seed=1572362164
Configuration cubic_side_length=$length particle_type=lj:/feasst/particle/lj.txt
Potential Model=LennardJones VisitModel=VisitModelCell
Potential VisitModel=LongRangeCorrections
Checkpoint checkpoint_file=checkpoint.fst
ThermoParams beta=$beta chemical_potential=-1
Metropolis
TrialTranslate weight=1 tunable_param=2
Tune
CheckEnergy trials_per_update=$tpc decimal_places=8
Log trials_per_write=$tpc output_file=lj_eq.csv
Run until_num_particles=$num_particles particle_type=lj Trial=TrialAdd weight=2
Run num_trials=1e5
Remove name=Tune,Log
WriteCheckpoint
EOF
