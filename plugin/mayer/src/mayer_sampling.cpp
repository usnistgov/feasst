#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h" // factorial
#include "math/include/constants.h"
#include "math/include/random.h"
#include "mayer/include/mayer_sampling.h"

namespace feasst {

MayerSampling::MayerSampling(argtype * args) : Criteria(args) {
  class_name_ = "MayerSampling";
  intra_pot_ = integer("intra_potential", args, -1);

  // HWH depreciate
  // Support depreciation warning for old argument name
  if (used("num_attempts_per_iteration", *args)) {
    WARN("Metropolis argument num_attempts_per_iteration is depreciated. " <<
         "Use num_trials_per_iteration instead.");
    ASSERT(!used("num_trials_per_iteration", *args),
      "Both num_trials_per_iteration and num_attempts_per_iteration");
    num_trials_per_iteration_ =
      integer("num_attempts_per_iteration", args);
  }
  num_trials_per_iteration_ =
    integer("num_trials_per_iteration", args, 1e9);
  const int num_beta_taylor = integer("num_beta_taylor", args, 0);
  beta_taylor_.resize(num_beta_taylor);
  for (int ibt = 0; ibt < num_beta_taylor; ++ibt) {
    beta_taylor_[ibt] = *MakeAccumulator({{"num_moments", "2"},
                                          {"max_block_operations", "0"}});
  }
}
MayerSampling::MayerSampling(argtype args) : MayerSampling(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void MayerSampling::precompute(System * system) {
  system->remove_opt_overlap();
}

bool MayerSampling::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  check_num_iterations_(num_trials_per_iteration_);
  ASSERT(system.num_configurations() == 1, "assumes 1 config");
  double energy_new = acceptance->energy_new();
  const double beta = system.thermo_params().beta();
  TRACE("*** MayerSampling ***");
  if (intra_pot_ != -1) {
    TRACE("intra_pot " << intra_pot_);
    const double energy_intra = acceptance->energy_profile_new()[intra_pot_];
    TRACE("energy_intra " << energy_intra);
    if (random->uniform() < std::exp(-beta*energy_intra)) {
      energy_new -= energy_intra;
    } else {
      was_accepted_ = false;
      TRACE("rejected at intra_potential step");
      return was_accepted_;
    }
  }
  const double f12 = std::exp(-beta*energy_new) - 1.;
  TRACE("energy new " << energy_new);
  TRACE("f12 " << f12);
  TRACE("f12old " << f12old_);
  TRACE("acceptance " << std::abs(f12)/std::abs(f12old_));

  if (!acceptance->reject() &&
      (random->uniform() < std::abs(f12)/std::abs(f12old_))) {
    ASSERT(energy_new != 0, "error");
    set_current_energy(acceptance->energy_new());
    set_current_energy_profile(acceptance->energy_profile_new());
    f12old_ = f12;
    was_accepted_ = true;
    TRACE("computing ref");
    f12ref_ = std::exp(-beta*acceptance->energy_ref()) - 1.;
    TRACE("f12ref " << f12ref_);
  } else {
    was_accepted_ = false;
  }
  if (f12old_ < 0) {
    mayer_.accumulate(-1.);
  } else {
    mayer_.accumulate(1.);
  }
  const double mayer_pi = std::abs(f12old_);
  mayer_ref_.accumulate(f12ref_/mayer_pi);
  double energy_old = current_energy();
  if (intra_pot_ != -1) {
    energy_old -= current_energy_profile()[intra_pot_];
  }
  double unebu = std::exp(-beta*energy_old)/mayer_pi;
  for (int ibd = 0; ibd < static_cast<int>(beta_taylor_.size()); ++ibd) {
    unebu *= -energy_old;
    beta_taylor_[ibd].accumulate(unebu/factorial(ibd+1));
  }
  TRACE("is accepted? " << was_accepted_);
  return was_accepted_;
}

class MapMayerSampling {
 public:
  MapMayerSampling() {
    MayerSampling().deserialize_map()["MayerSampling"] = MakeMayerSampling();
  }
};

static MapMayerSampling mapper_ = MapMayerSampling();

void MayerSampling::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(3252, ostr);
  feasst_serialize(f12old_, ostr);
  feasst_serialize(f12ref_, ostr);
  feasst_serialize(num_trials_per_iteration_, ostr);
  feasst_serialize_fstobj(mayer_, ostr);
  feasst_serialize_fstobj(mayer_ref_, ostr);
  feasst_serialize(intra_pot_, ostr);
  feasst_serialize_fstobj(beta_taylor_, ostr);
}

MayerSampling::MayerSampling(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3251 && version <= 3252, "unrecognized verison: " << version);
  feasst_deserialize(&f12old_, istr);
  feasst_deserialize(&f12ref_, istr);
  feasst_deserialize(&num_trials_per_iteration_, istr);
  feasst_deserialize_fstobj(&mayer_, istr);
  feasst_deserialize_fstobj(&mayer_ref_, istr);
  feasst_deserialize(&intra_pot_, istr);
  if (version >= 3252) {
    feasst_deserialize_fstobj(&beta_taylor_, istr);
  }
}

double MayerSampling::second_virial_ratio() const {
  return mayer_.average()/mayer_ref_.average();
}

//MayerSampling::MayerSampling(const Criteria& criteria) {
//  std::stringstream ss;
//  criteria.serialize(ss);
//  *this = MayerSampling(ss);
//}

double MayerSampling::second_virial_ratio_block_stdev() const {
  return std::sqrt(
    pow(mayer_.block_stdev()/mayer_ref_.average(), 2) +
    pow(second_virial_ratio()*mayer_ref_.block_stdev()/mayer_ref_.average(), 2));
}

std::string MayerSampling::write() const {
  std::stringstream ss;
  ss << "{\"second_virial_ratio\": " << second_virial_ratio() << ", "
     << "\"second_virial_ratio_block_stdev\": " << second_virial_ratio_block_stdev();
  if (num_beta_taylor() > 0) {
    ss << ",\"beta_taylor\": [";
    for (int ibt = 0; ibt < num_beta_taylor() + 1; ++ibt) {
      ss << beta_taylor(ibt) << ",";
    }
    ss << "],";
  }
  ss << "}" << std::endl;
  ss << "mayer" << std::endl << mayer().str() << std::endl;
  ss << "mayer_ref" << std::endl << mayer_ref().str() << std::endl;
  return ss.str();
}

double MayerSampling::beta_taylor(const int deriv) const {
  if (deriv == 0) {
    return mayer().average()/mayer_ref().average();
  } else {
    return beta_taylor_[deriv -1].average()/mayer_ref().average();
  }
}

}  // namespace feasst
