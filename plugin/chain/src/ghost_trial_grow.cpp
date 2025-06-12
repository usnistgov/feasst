#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/random_mt19937.h" // HWH remove this
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_factory.h"
#include "chain/include/trial_grow.h"
#include "chain/include/ghost_trial_grow.h"

namespace feasst {

FEASST_MAPPER(GhostTrialGrow,);

GhostTrialGrow::GhostTrialGrow() : Modify() {} // only use for deserialize_map.
GhostTrialGrow::~GhostTrialGrow() {}

void GhostTrialGrow::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2948, ostr);
  feasst_serialize(grow_, ostr);
  feasst_serialize_fstobj(criteria_, ostr);
  feasst_serialize_fstobj(metropolis_prob_, ostr);
}

GhostTrialGrow::GhostTrialGrow(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2948 == version, version);
  feasst_deserialize(grow_, istr);
  feasst_deserialize_fstobj(&criteria_, istr);
  feasst_deserialize_fstobj(&metropolis_prob_, istr);
}

GhostTrialGrow::GhostTrialGrow(argtype * args) : Modify(args) {
  const std::string grow_file = str("grow_file", args);
  auto grow = MakeTrialGrowFile({{"grow_file", grow_file}});
  const int num_trials = static_cast<int>(grow->trials().size());
  ASSERT(num_trials == 1, "GhostTrialGrow needs just one trial, but there are "
    << num_trials << " trials.");
  grow_ = std::make_unique<TrialFactory>();
  grow_->add(grow->trials()[0]);
}
GhostTrialGrow::GhostTrialGrow(argtype args) : GhostTrialGrow(&args) {
  feasst_check_all_used(args);
}

void GhostTrialGrow::initialize(MonteCarlo * mc) {
  Modify::initialize(mc);
  System * system = mc->get_system();
  Criteria * criteria = mc->get_criteria();
  printer(header(*mc), output_file(mc->criteria()));
  // HWH double check this initialization
  const int conf = configuration_index();
  criteria_.set_current_energy(criteria->current_energy(), conf);
  criteria_.set_current_energy_profile(system->stored_energy_profile(conf), conf);
  criteria_.precompute(system);
  grow_->get_trial(0)->precompute(&criteria_, system);
}

void GhostTrialGrow::update(MonteCarlo * mc) {
  System * system = mc->get_system();
  Random * random = mc->get_random();
  const int trial_index = 0;
  grow_->attempt(&criteria_, system, trial_index, random);
  const Acceptance& acc = grow_->trial(0).accept();
  //const int conf = configuration_index();
  //const double beta = system->thermo_params().beta();
  //INFO(acc.ln_metropolis_prob() << " "
  //  << acc.energy_new(conf) << " "
  //  << acc.energy_old(conf));
  //const double delta_energy = acc.energy_new(conf) - acc.energy_old(conf);
  if (acc.reject()) {
    metropolis_prob_.accumulate(0.);
  } else {
    metropolis_prob_.accumulate(std::exp(acc.ln_metropolis_prob()));
  }
  grow_->revert(trial_index, false, acc.reject(), system, &criteria_);
}

std::string GhostTrialGrow::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << metropolis_prob_.status_header() << std::endl;
  return ss.str();
}

std::string GhostTrialGrow::write(MonteCarlo * mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*mc);
  }
  ss << metropolis_prob_.status() << std::endl;
  return ss.str();
}

}  // namespace feasst
