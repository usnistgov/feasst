#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h" // factorial
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_stage.h"
#include "confinement/include/henry_coefficient.h"

namespace feasst {

class MapHenryCoefficient {
 public:
  MapHenryCoefficient() {
    auto obj = MakeHenryCoefficient();
    obj->deserialize_map()["HenryCoefficient"] = obj;
  }
};

static MapHenryCoefficient mapper_ = MapHenryCoefficient();

HenryCoefficient::HenryCoefficient(argtype * args) : Analyze(args) {
  ASSERT(trials_per_update() == 1, "should update every step");
  const int num_beta_taylor = integer("num_beta_taylor", args, 0);
  beta_taylor_.resize(num_beta_taylor);
  for (int ibt = 0; ibt < num_beta_taylor; ++ibt) {
    beta_taylor_[ibt] = *MakeAccumulator({{"num_moments", "2"},
                                          {"max_block_operations", "0"}});
  }
}
HenryCoefficient::HenryCoefficient(argtype args) : HenryCoefficient(&args) {
  feasst_check_all_used(args);
}

void HenryCoefficient::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  ASSERT(criteria->class_name() == "AlwaysReject",
    "HenryCoefficient requires AlwaysReject");
  ASSERT(trial_factory->num() == 1,
    "HenryCoefficient requires only one Trial");
  ASSERT(trial_factory->trial(0).description() == "TrialAdd",
    "HenryCoefficient requires TrialAdd. " <<
    "Found: " << trial_factory->trial(0).description());
  ASSERT(trial_factory->trial(0).stage(0).is_new_only(),
    "HenryCoefficient requires new_only.");
  printer(header(*criteria, *system, *trial_factory),
          output_file(*criteria));
}

std::string HenryCoefficient::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void HenryCoefficient::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const double en = trial_factory.trial(0).accept().energy_new();
  DEBUG("en: " << en);
  const double beta = system.thermo_params().beta();
  get_accumulator()->accumulate(std::exp(-beta*en));

  // moments for extrapolation
  double unebu = std::exp(-beta*en);
  for (int ibd = 0; ibd < static_cast<int>(beta_taylor_.size()); ++ibd) {
    unebu *= -en;
    beta_taylor_[ibd].accumulate(unebu/factorial(ibd+1));
  }
}

std::string HenryCoefficient::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << accumulator().status();
  if (num_beta_taylor() > 0) {
    ss << ",\"beta_taylor\": [";
    ss << accumulator().average() << ",";
    for (int ibt = 1; ibt < num_beta_taylor() + 1; ++ibt) {
      ss << beta_taylor_[ibt - 1].average() << ",";
    }
    ss << "],";
  }
  ss << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void HenryCoefficient::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(9493, ostr);
  feasst_serialize_fstobj(beta_taylor_, ostr);
}

HenryCoefficient::HenryCoefficient(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 9492 && version <= 9493, "mismatch version:" << version);
  if (version >= 9493) {
    feasst_deserialize_fstobj(&beta_taylor_, istr);
  }
}

}  // namespace feasst
