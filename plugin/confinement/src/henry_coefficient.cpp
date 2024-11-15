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
#include "monte_carlo/include/monte_carlo.h"
#include "confinement/include/henry_coefficient.h"

namespace feasst {

FEASST_MAPPER(HenryCoefficient,);

HenryCoefficient::HenryCoefficient(argtype * args) : Analyze(args) {
  ASSERT(trials_per_update() == 1, "should update every step");
  const int num_beta_taylor = integer("num_beta_taylor", args, 0);
  beta_taylor_.resize(num_beta_taylor);
  beta_taylor2_.resize(2*num_beta_taylor);
  for (int ibt = 0; ibt < num_beta_taylor; ++ibt) {
    beta_taylor_[ibt] = *MakeAccumulator({{"num_moments", "2"},
                                          {"max_block_operations", "0"}});
  }
  for (int ibt = 0; ibt < 2*num_beta_taylor; ++ibt) {
    beta_taylor2_[ibt] = *MakeAccumulator({{"num_moments", "2"},
                                          {"max_block_operations", "0"}});
  }
  write_precision_ = dble("write_precision", args, 8);
}
HenryCoefficient::HenryCoefficient(argtype args) : HenryCoefficient(&args) {
  feasst_check_all_used(args);
}

void HenryCoefficient::initialize(MonteCarlo * mc) {
  const Criteria& criteria = mc->criteria();
  const TrialFactory& tfac = mc->trial_factory();
  ASSERT(criteria.class_name() == "AlwaysReject",
    "HenryCoefficient requires AlwaysReject");
  ASSERT(tfac.num() == 1,
    "HenryCoefficient requires only one Trial");
  ASSERT(tfac.trial(0).description() == "TrialAdd",
    "HenryCoefficient requires TrialAdd. " <<
    "Found: " << tfac.trial(0).description());
  ASSERT(tfac.trial(0).stage(0).is_new_only(),
    "HenryCoefficient requires new_only.");
  printer(header(*mc), output_file(mc->criteria()));
}

std::string HenryCoefficient::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void HenryCoefficient::update(const MonteCarlo& mc) {
  const Acceptance& accept = mc.trial_factory().trial(0).accept();
  double en_new = accept.energy_new();
  DEBUG("en new: " << en_new);
  double en_cur = mc.criteria().current_energy();
  DEBUG("en cur: " << en_cur);
  const double delta_en = en_new - en_cur;
  DEBUG("delta_en " << delta_en);
  const double beta = mc.system().thermo_params().beta();
  DEBUG("beta " << beta);
  get_accumulator()->accumulate(std::exp(-beta*delta_en));
  DEBUG("acc " << accumulator().str());

  // moments for extrapolation including covariance terms
  double unebu = std::exp(-beta*delta_en);
  double unebu2 = std::exp(-2.0*beta*delta_en);
  DEBUG("unebu " << unebu);
  DEBUG("unebu2 " << unebu2);
  for (int ibd = 0; ibd < static_cast<int>(beta_taylor_.size()); ++ibd) {
    unebu *= -delta_en;
    DEBUG("ibd " << ibd << " unebu " << unebu);
    beta_taylor_[ibd].accumulate(unebu/factorial(ibd+1));
    unebu2 *= -delta_en;
    DEBUG("2*ibd " << 2*ibd << " unebu2 " << unebu2);
    beta_taylor2_[2*ibd].accumulate(unebu2/factorial(2*ibd+1)/factorial(2*ibd+1));
    unebu2 *= -delta_en;
    DEBUG("2*ibd+1 " << 2*ibd+1 << " unebu2 " << unebu2);
    beta_taylor2_[2*ibd+1].accumulate(unebu2/factorial(2*ibd+2)/factorial(2*ibd+2));
  }
}

std::string HenryCoefficient::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (num_beta_taylor() > 0) {
    ss << "#{";
    ss << "\"beta\":" << std::setprecision(write_precision_) << mc.system().thermo_params().beta() << ",";
    ss << "\"num_trials\":" << accumulator().moment(0) << ",";
    ss << "\"beta_taylor\": [";
    ss << std::setprecision(write_precision_) << accumulator().average() << ",";
    for (int ibt = 1; ibt < num_beta_taylor() + 1; ++ibt) {
      ss << std::setprecision(write_precision_)
         << beta_taylor_[ibt - 1].average();
      if (ibt < num_beta_taylor()) ss << ",";
    }
    ss << "], \"beta_taylor2\": [";
    ss << std::setprecision(write_precision_) << accumulator().moment(2)/accumulator().moment(0) << ",";
    for (int ibt = 1; ibt < 2*num_beta_taylor() + 1; ++ibt) {
      ss << std::setprecision(write_precision_)
         << beta_taylor2_[ibt - 1].average();
      if (ibt < 2*num_beta_taylor()) ss << ",";
    }
    ss << "]}";
  }
  ss << std::endl;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator().status();
  DEBUG(ss.str());
  return ss.str();
}

void HenryCoefficient::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(9494, ostr);
  feasst_serialize_fstobj(beta_taylor_, ostr);
  feasst_serialize_fstobj(beta_taylor2_, ostr);
  feasst_serialize(write_precision_, ostr);
}

HenryCoefficient::HenryCoefficient(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 9492 && version <= 9494, "mismatch version:" << version);
  if (version >= 9493) {
    feasst_deserialize_fstobj(&beta_taylor_, istr);
  }
  if (version >= 9494) {
    feasst_deserialize_fstobj(&beta_taylor2_, istr);
    feasst_deserialize(&write_precision_, istr);
  }
}

}  // namespace feasst
