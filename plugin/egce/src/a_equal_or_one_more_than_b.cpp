#include "utils/include/serialize.h"
#include "egce/include/a_equal_or_one_more_than_b.h"

namespace feasst {

bool AEqualOrOneMoreThanB::is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const {
  const int na = system->configuration().num_particles_of_type(0);
  const int nb = system->configuration().num_particles_of_type(1);
  bool allowed = false;
  const int shift = acceptance.macrostate_shift();
  if (shift == 0) {
    if (na == nb || na == nb + 1) {
      allowed = true;
    }
  } else if (std::abs(shift) == 1) {
    if (acceptance.macrostate_shift_type() == 0) {
      if (na == nb + 1) {
        allowed = true;
      }
    } else if (acceptance.macrostate_shift_type() == 1) {
      if (na == nb) {
        allowed = true;
      }
    } else {
      FATAL("unrecognized shift type:" << acceptance.macrostate_shift_type());
    }
  } else {
    FATAL("unrecognized shift: " << shift);
  }
  DEBUG("na " << na << " nb " << nb << " allowed " << allowed << " shift "
    << shift << " type " << acceptance.macrostate_shift_type());
  return allowed;
}

class MapAEqualOrOneMoreThanB {
 public:
  MapAEqualOrOneMoreThanB() {
    auto obj = MakeAEqualOrOneMoreThanB();
    obj->deserialize_map()["AEqualOrOneMoreThanB"] = obj;
  }
};

static MapAEqualOrOneMoreThanB mapper_ = MapAEqualOrOneMoreThanB();

std::shared_ptr<Constraint> AEqualOrOneMoreThanB::create(std::istream& istr) const {
  return std::make_shared<AEqualOrOneMoreThanB>(istr);
}

AEqualOrOneMoreThanB::AEqualOrOneMoreThanB(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "AEqualOrOneMoreThanB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2492 == version, "mismatch version: " << version);
}

void AEqualOrOneMoreThanB::serialize_a_equal_or_one_more_than_b_(std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2492, ostr);
}

void AEqualOrOneMoreThanB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_a_equal_or_one_more_than_b_(ostr);
}

}  // namespace feasst
