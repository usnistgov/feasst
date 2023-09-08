#include <cmath>
#include "utils/include/serialize.h"
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
  mayer_ref_.accumulate(f12ref_/std::abs(f12old_));
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
  feasst_serialize_version(3251, ostr);
  feasst_serialize(f12old_, ostr);
  feasst_serialize(f12ref_, ostr);
  feasst_serialize(num_trials_per_iteration_, ostr);
  feasst_serialize_fstobj(mayer_, ostr);
  feasst_serialize_fstobj(mayer_ref_, ostr);
  feasst_serialize(intra_pot_, ostr);
}

MayerSampling::MayerSampling(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3251, "unrecognized verison: " << version);
  feasst_deserialize(&f12old_, istr);
  feasst_deserialize(&f12ref_, istr);
  feasst_deserialize(&num_trials_per_iteration_, istr);
  feasst_deserialize_fstobj(&mayer_, istr);
  feasst_deserialize_fstobj(&mayer_ref_, istr);
  feasst_deserialize(&intra_pot_, istr);
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
     << "\"second_virial_ratio_block_stdev\": " << second_virial_ratio_block_stdev() << "}" << std::endl;
  ss << "mayer" << std::endl << mayer().str() << std::endl;
  ss << "mayer_ref" << std::endl << mayer_ref().str() << std::endl;
  return ss.str();
}

}  // namespace feasst
