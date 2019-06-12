#include <iostream>
#include "feasst.h"

int main() {
  feasst::seed_random_by_date();
  feasst::MonteCarlo mc;
  mc.add(feasst::Configuration({
    {"cubic_box_length", "8"},
    {"particle_type", "../../../../forcefield/data.lj"},
  }));
  mc.add(feasst::Potential(feasst::MakeModelLJ()));
  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  mc.add(feasst::MakeCriteriaMetropolis({
    {"beta", "1.2"},
    {"chemical_potential", "1."},
  }));
  mc.add(feasst::MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "2."},
  }));
  mc.seek_num_particles(50);
  const int steps_per = 1e4;
  mc.add(feasst::MakeLog({{"steps_per", feasst::str(steps_per)}}));
  mc.add(feasst::MakeMovie(
   {{"steps_per", feasst::str(steps_per)},
    {"file_name", "movie.xyz"}}));
  mc.add(feasst::MakeCheckEnergy(
   {{"steps_per", feasst::str(steps_per)},
    {"tolerance", "1e-10"}}));
  mc.add(feasst::MakeTuner({{"steps_per", feasst::str(steps_per)}}));
  mc.attempt(1e6);
}
