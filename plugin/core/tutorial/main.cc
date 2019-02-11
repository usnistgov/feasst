#include <iostream>
#include "feasst.h"

int main() {
  feasst::seed_random_by_date();
  feasst::MonteCarlo mc;

  { // add system to mc
    feasst::System sys;
    { // add configuration to system
      feasst::Configuration config;
      config.set_domain(feasst::Domain().set_cubic(8));
      config.add_particle_type("../../../../forcefield/data.lj");
      sys.add_configuration(config);
    }

    { // add potentials to system
      feasst::Potential potential;
      potential.set_model(std::make_shared<feasst::ModelLJ>());
      potential.set_visit_model(std::make_shared<feasst::VisitModel>());
      feasst::Potentials potentials;
      potentials.add_potential(potential);
      potential.set_model(std::make_shared<feasst::ModelLRC>());
      potentials.add_potential(potential);
      sys.set_full(potentials);
    }
    mc.set_system(sys);
  }

  { // add criteria to mc
    auto criteria = std::make_shared<feasst::CriteriaMetropolis>();
    criteria->set_beta(1.2);
    criteria->add_activity(1);
    mc.set_criteria(criteria);
  }


  { // add translate to mc
    auto translate = std::make_shared<feasst::TrialTranslate>();
    translate->set_weight(1);
    mc.add_trial(translate);
  }

  mc.seek_num_particles(50);

  { // add analyze to mc
    auto log = std::make_shared<feasst::Log>();
    log->set_steps_per_write(1e4);
    mc.add_analyze(log);
  }

  mc.attempt(1e6);
}
