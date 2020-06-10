#include "feasst.h"

int main() {
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", "time"}}));
  mc.add(feasst::Configuration(feasst::MakeDomain({{"cubic_box_length", "8"}}),
    {{"particle_type", feasst::install_dir() + "/forcefield/data.lj"}}));
  mc.add(feasst::Potential(feasst::MakeLennardJones()));
  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  mc.add(feasst::MakeMetropolis({{"beta", "1.2"}}));
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
  const int steps_per = 1e3;
  mc.add(feasst::MakeCheckEnergyAndTune(
   {{"steps_per", feasst::str(steps_per)}, {"tolerance", "1e-8"}}));
  feasst::SeekNumParticles(50)
    .with_metropolis({{"beta", "0.1"}, {"chemical_potential", "10"}})
    .with_trial_add()
    .run(&mc);
  mc.add(feasst::MakeLogAndMovie(
   {{"steps_per", feasst::str(steps_per)}, {"file_name", "lj"}}));
  mc.attempt(1e5);
}
