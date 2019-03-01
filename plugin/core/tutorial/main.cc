#include <iostream>
#include "feasst.h"

feasst::Configuration configuration() {
  feasst::Configuration config;
  config.set_domain(feasst::Domain().set_cubic(8));
  config.add_particle_type("../../../../forcefield/data.lj");
  return config;
}

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
  system.add(configuration());
  system.add(lj());
  system.add(lrc());
  return system;
}

std::shared_ptr<feasst::Criteria> criteria(feasst::System * system) {
  auto criteria = std::make_shared<feasst::CriteriaMetropolis>();
  criteria->set_beta(1.2);
  criteria->add_activity(1);
  criteria->set_running_energy(system->energy());
  return criteria;
}

std::shared_ptr<feasst::Trial> translate(const feasst::Domain& domain) {
  auto trial = std::make_shared<feasst::TrialTranslate>();
  trial->set_weight(1);
  trial->set_max_move(2.);
  return trial;
}

int main() {
  feasst::seed_random_by_date();
  feasst::MonteCarlo mc;
  mc.set(system());
  mc.set(criteria());
  mc.add(translate(mc.system().configuration().domain()));
  mc.seek_num_particles(50);
  const int num_periodic = 1e4;
  mc.add(feasst::MakeLog({{"steps_per", feasst::str(num_periodic)}}));
  mc.add(feasst::MakeMovie(
   {{"steps_per", feasst::str(num_periodic)},
    {"file_name", "tmp/lj50movie.xyz"}}));
  mc.add(feasst::MakeEnergyCheck(
   {{"steps_per", feasst::str(num_periodic)},
    {"tolerance", "1e-10"}}));
  mc.add(feasst::MakeTuner({{"steps_per", feasst::str(num_periodic)}}));
  mc.attempt(1e6);
}
