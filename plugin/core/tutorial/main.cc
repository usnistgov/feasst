#include <iostream>
#include "feasst.h"

feasst::Potential lj() {
  feasst::Potential potential;
  potential.set_model(std::make_shared<feasst::ModelLJ>());
  potential.set_visit_model(std::make_shared<feasst::VisitModel>());
  return potential;
}

feasst::Potential lrc() {
  feasst::Potential potential;
  potential.set_visit_model(std::make_shared<feasst::LongRangeCorrections>());
  return potential;
}

feasst::System system() {
  feasst::System system;
  system.add(feasst::Configuration({
    {"cubic_box_length", "8"},
    {"particle_type", "../../../../forcefield/data.lj"},
  }));
  system.add(lj());
  system.add(lrc());
  return system;
}

int main() {
  feasst::seed_random_by_date();
  feasst::MonteCarlo mc;
  mc.set(system());
  mc.set(feasst::MakeCriteriaMetropolis({
    {"beta", "1.2"},
    {"add_activity", "1."},
  }));
  mc.add(feasst::MakeTrialTranslate({
    {"weight", "1."},
    {"max_move", "2."},
  }));
  mc.seek_num_particles(50);
  const int steps_per = 1e4;
  mc.add(feasst::MakeLog({{"steps_per", feasst::str(steps_per)}}));
  mc.add(feasst::MakeMovie(
   {{"steps_per", feasst::str(steps_per)},
    {"file_name", "movie.xyz"}}));
  mc.add(feasst::MakeEnergyCheck(
   {{"steps_per", feasst::str(steps_per)},
    {"tolerance", "1e-10"}}));
  mc.add(feasst::MakeTuner({{"steps_per", feasst::str(steps_per)}}));
  mc.attempt(1e6);
}
