import sys
import subprocess
import unittest
import feasst as fst
import pyfeasst

class TestMonteCarlo1LJNVT(unittest.TestCase):
    """Test a canonical ensemble Lennard Jones Monte Carlo simulation"""
    def test_srsw_alt(self, num_particles=500, density=0.001, steps_per=1e5, beta=1./0.9,
                      fstprt=fst.install_dir() + "/forcefield/lj.fstprt",
                      num_equil=1e7, num_prod=1e7):
        """Compare with the reported average energy from the NIST SRSW.
        https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
        https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website

        num_particles -- number of LJ particles
        density -- number density
        steps_per -- steps between each Anaylze/Modify
        """
        params = {"box_length": (num_particles/density)**(1./3.)}
        params = dict(locals(), **params)
        with open('tutorial_1_lj_nvt.txt', 'w') as fsttxt:
            fsttxt.write("""
Checkpoint file_name checkpoint.fst
RandomMT19937 seed time
Configuration cubic_box_length {box_length} particle_type {fstprt}
Potential Model LennardJones
Potential VisitModel LongRangeCorrections
ThermoParams beta 0.1 chemical_potential 10
Metropolis
TrialTranslate tunable_param 2. tunable_target_acceptance 0.2
TrialAdd particle_type 0
Run until_num_particles {num_particles}
RemoveTrial name TrialAdd
ThermoParams beta {beta}
Tune steps_per {steps_per}
CheckEnergy steps_per {steps_per} tolerance 1e-8

# equilibrate
Run num_attempts {num_equil}
RemoveModify name Tune

# production analysis and output
Log steps_per {steps_per} file_name lj.txt
Energy steps_per_write {steps_per} file_name en.txt
# Run num_attempts {num_prod} # optionally run production here instead of python
WriteCheckpoint
""".format(**params))
        import pyfeasst
        syscode = subprocess.call(fst.install_dir() + "/build/bin/fst < tutorial_1_lj_nvt.txt >> launch.log", shell=True, executable='/bin/bash')
        if syscode > 0: sys.exit(1)
        mc = fst.MonteCarlo().deserialize(pyfeasst.read_checkpoint('checkpoint.fst'))

        # run production trials with an alternative running average of the energy
        # this demonstrates custom on-the-fly analysis in python scripts.
        energy_alt = fst.Accumulator()
        for _ in range(int(params["num_prod"])):
            mc.attempt(1)
            energy_alt.accumulate(mc.criteria().current_energy())
        energy = fst.SeekAnalyze().reference("Energy", mc).accumulator()

        # test that the two methods for computing the energy give the same result
        self.assertAlmostEqual(energy.average(), energy_alt.average(), delta=1e-6)

        # test the average against the NIST SRSW
        stdev = (energy.block_stdev()**2 + (1.89E-05)**2)**(1./2.)
        self.assertAlmostEqual(-9.9165E-03*num_particles, energy.average(),
                               delta=2.576*stdev)

unittest.main(argv=[''], verbosity=2, exit=False)
