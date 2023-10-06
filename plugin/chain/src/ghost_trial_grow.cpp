#include <cmath>
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/random_mt19937.h" // HWH remove this
#include "configuration/include/domain.h"
#include "system/include/ideal_gas.h"
#include "chain/include/trial_grow.h"
#include "chain/include/ghost_trial_grow.h"

namespace feasst {

class MapGhostTrialGrow {
 public:
  MapGhostTrialGrow() {
    auto obj = std::make_shared<GhostTrialGrow>();
    //auto obj = MakeGhostTrialGrow({{"trial_grow_file", install_dir()+"/plugin/chain/test/data/dimer_grow_file.txt"}});
    obj->deserialize_map()["GhostTrialGrow"] = obj;
  }
};

static MapGhostTrialGrow mapper_inner_ = MapGhostTrialGrow();

void GhostTrialGrow::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2948, ostr);
  feasst_serialize_fstobj(metropolis_prob_, ostr);
  feasst_serialize_fstobj(grow_, ostr);
}

GhostTrialGrow::GhostTrialGrow(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2948 == version, version);
  feasst_deserialize_fstobj(&metropolis_prob_, istr);
  feasst_deserialize_fstobj(&grow_, istr);
}

GhostTrialGrow::GhostTrialGrow(argtype * args) : Modify(args) {
  const std::string trial_grow_file = str("trial_grow_file", args);
  auto grow = MakeTrialGrowFile({{"file_name", trial_grow_file}});
  const int num_trials = static_cast<int>(grow->trials().size());
  ASSERT(num_trials == 1, "GhostTrialGrow needs just one trial, but there are "
    << num_trials << " trials.");
  grow_.add(grow->trials()[0]);
}
GhostTrialGrow::GhostTrialGrow(argtype args) : GhostTrialGrow(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void GhostTrialGrow::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
  const int conf = configuration_index();
  criteria_.set_current_energy(criteria->current_energy(), conf);
  criteria_.set_current_energy_profile(system->stored_energy_profile(conf), conf);
  criteria_.precompute(system);
  grow_.get_trial(0)->precompute(&criteria_, system);
}

void GhostTrialGrow::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  const int trial_index = 0;
  grow_.attempt(&criteria_, system, trial_index, random);
  const Acceptance& acc = grow_.trial(0).accept();
  //const int conf = configuration_index();
  //const double beta = system->thermo_params().beta();
  //INFO(acc.ln_metropolis_prob() << " "
  //  << acc.energy_new(conf) << " "
  //  << acc.energy_old(conf));
  //const double delta_energy = acc.energy_new(conf) - acc.energy_old(conf);
  metropolis_prob_.accumulate(std::exp(acc.ln_metropolis_prob()));
  grow_.revert(trial_index, false, acc.reject(), system, &criteria_);
}

std::string GhostTrialGrow::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << metropolis_prob_.status_header() << std::endl;
  return ss.str();
}

std::string GhostTrialGrow::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*criteria, *system, *trial_factory);
  }
  ss << metropolis_prob_.status() << std::endl;
  return ss.str();
}

}  // namespace feasst
