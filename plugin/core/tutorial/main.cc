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

std::shared_ptr<feasst::Criteria> criteria() {
  auto criteria = std::make_shared<feasst::CriteriaMetropolis>();
  criteria->set_beta(1.2);
  criteria->add_activity(1);
  return criteria;
}

std::shared_ptr<feasst::Trial> translate(const feasst::Domain& domain) {
  auto trial = std::make_shared<feasst::TrialTranslate>();
  trial->set_weight(1);
  trial->set_max_move(2.);
  trial->set_max_move_bounds(domain);
  return trial;
}

std::shared_ptr<feasst::Log> log() {
  auto log = std::make_shared<feasst::Log>();
  log->set_steps_per_write(1e4);
  return log;
}

std::shared_ptr<feasst::Movie> movie() {
  auto movie = std::make_shared<feasst::Movie>();
  movie->set_steps_per_write(1e4);
  movie->set_file_name("movie.xyz");
  return movie;
}

std::shared_ptr<feasst::EnergyCheck> check() {
  auto check = std::make_shared<feasst::EnergyCheck>();
  check->set_steps_per_update(1e4);
  check->set_tolerance(1e-10);
  return check;
}

std::shared_ptr<feasst::Tuner> tune() {
  auto tune = std::make_shared<feasst::Tuner>();
  tune->set_steps_per_update(1e4);
  return tune;
}

int main() {
  feasst::seed_random_by_date();
  feasst::MonteCarlo mc;
  mc.set(system());
  mc.set(criteria());
  mc.add(translate(mc.system().configuration().domain()));
  mc.seek_num_particles(50);
  mc.add(log());
  mc.add(movie());
  mc.add(check());
  mc.add(tune());
  mc.attempt(1e6);
}
