#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

void MonteCarlo::add(std::shared_ptr<Analyze> analyze) {
  ASSERT(criteria_set_, "set Criteria before Analyze");
  DEBUG("class name? " << analyze->class_name());
  if (analyze->is_multistate() and analyze->class_name() != "AnalyzeFactory") {
    auto multi = MakeAnalyzeFactory({{"multistate", "true"}});
    DEBUG("making multi " << multi->is_multistate());
    for (int state = 0; state < criteria_->num_states(); ++state) {
      DEBUG("state: " << state);
      std::shared_ptr<Analyze> an = deep_copy_derived(analyze);
      { std::stringstream ss;
        an->serialize(ss);
        DEBUG(ss.str());
      }
      an->set_state(state);
      // an->initialize(criteria_, system_, trial_factory_);
      multi->add(an);
    }
    analyze = multi;
  }
  analyze->initialize(criteria_, system_, trial_factory_);
  DEBUG("mults " << analyze->is_multistate() << " class name? " << analyze->class_name());
  analyze_factory_.add(analyze);
}

}  // namespace feasst
