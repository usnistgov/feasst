#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "mayer/include/mayer_sampling.h"

namespace feasst {

MayerSampling::MayerSampling() : Criteria() {
  class_name_ = "MayerSampling";
}

bool MayerSampling::is_accepted(const Acceptance& acceptance,
    const System& system,
    Random * random) {
  const double energy_new = acceptance.energy_new();
  const double beta = system.thermo_params().beta();
  const double f12 = std::exp(-beta*energy_new) - 1.;
  bool is_accepted;
  DEBUG("*** MayerSampling ***");
  DEBUG("energy new " << energy_new);
  DEBUG("f12 " << f12);
  DEBUG("f12old " << f12old_);
  DEBUG("acceptance " << std::abs(f12)/std::abs(f12old_));
  if (!acceptance.reject() &&
      (random->uniform() < std::abs(f12)/std::abs(f12old_))) {
    set_current_energy(energy_new);
    f12old_ = f12;
    is_accepted = true;
    DEBUG("computing ref");
    f12ref_ = std::exp(-beta*acceptance.energy_ref()) - 1.;
    DEBUG("f12ref " << f12ref_);
  } else {
    is_accepted = false;
  }
  if (f12old_ < 0) {
    mayer_.accumulate(-1.);
  } else {
    mayer_.accumulate(1.);
  }
  mayer_ref_.accumulate(f12ref_/std::abs(f12old_));
  DEBUG("is accepted? " << is_accepted);
  return is_accepted;
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
  feasst_serialize_fstobj(mayer_, ostr);
  feasst_serialize_fstobj(mayer_ref_, ostr);
}

MayerSampling::MayerSampling(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3251, "unrecognized verison: " << version);
  feasst_deserialize(&f12old_, istr);
  feasst_deserialize(&f12ref_, istr);
  feasst_deserialize_fstobj(&mayer_, istr);
  feasst_deserialize_fstobj(&mayer_ref_, istr);
}

double MayerSampling::second_virial_ratio() const {
  return mayer_.average()/mayer_ref_.average();
}

//MayerSampling::MayerSampling(const Criteria& criteria) {
//  std::stringstream ss;
//  criteria.serialize(ss);
//  *this = MayerSampling(ss);
//}

}  // namespace feasst
