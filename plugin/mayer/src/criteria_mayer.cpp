
#include "mayer/include/criteria_mayer.h"
#include "utils/include/utils_io.h"

namespace feasst {

bool CriteriaMayer::is_accepted(const Acceptance& acceptance,
    const System * system) {
  const double energy_new = acceptance.energy_new();
  const double f12 = exp(-beta()*energy_new) - 1.;
  bool is_accepted;
  DEBUG("energy new " << energy_new << " f12 " << f12);
  if (!acceptance.reject() and
      (random_.uniform() < std::abs(f12)/std::abs(f12old_)) ) {
    set_current_energy(energy_new);
    f12old_ = f12;
    is_accepted = true;
    DEBUG("computing ref");
    f12ref_ = exp(-beta()*acceptance.energy_ref()) - 1.;
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

class MapCriteriaMayer {
 public:
  MapCriteriaMayer() {
    CriteriaMayer().deserialize_map()["CriteriaMayer"] = MakeCriteriaMayer();
  }
};

static MapCriteriaMayer mapper_ = MapCriteriaMayer();

void CriteriaMayer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize(f12old_, ostr);
  feasst_serialize(f12ref_, ostr);
  feasst_serialize_fstobj(mayer_, ostr);
  feasst_serialize_fstobj(mayer_ref_, ostr);
}

CriteriaMayer::CriteriaMayer(std::istream& istr) : Criteria(istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&f12old_, istr);
  feasst_deserialize(&f12ref_, istr);
  feasst_deserialize_fstobj(&mayer_, istr);
  feasst_deserialize_fstobj(&mayer_ref_, istr);
}

}  // namespace feasst
