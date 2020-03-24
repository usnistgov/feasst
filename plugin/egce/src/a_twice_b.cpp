#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "egce/include/a_twice_b.h"

namespace feasst {

bool ATwiceB::is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const {
  int na = system->configuration().num_particles_of_type(0);
  int nb = system->configuration().num_particles_of_type(1);
  bool allowed = false;
  const int shift = acceptance.macrostate_shift();
  if (shift == -1) {
    const int type = acceptance.macrostate_shift_type();
    if (type == 0) {
      --na;
    } else if (type == 1) {
      --nb;
    } else {
      FATAL("unrecognized type: " << type);
    }
  } else {
    ASSERT(shift == 0, "unrecognized shift: " << shift);
  }
  if (std::abs(static_cast<double>(na)/2. - nb) <= 0.5 + NEAR_ZERO) {
    allowed = true;
  }
  DEBUG("na " << na << " nb " << nb << " allowed " << allowed << " shift "
    << shift << " type " << acceptance.macrostate_shift_type());
  return allowed;
}

class MapATwiceB {
 public:
  MapATwiceB() {
    auto obj = MakeATwiceB();
    obj->deserialize_map()["ATwiceB"] = obj;
  }
};

static MapATwiceB mapper_ = MapATwiceB();

std::shared_ptr<Constraint> ATwiceB::create(std::istream& istr) const {
  return std::make_shared<ATwiceB>(istr);
}

ATwiceB::ATwiceB(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "ATwiceB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2492 == version, "mismatch version: " << version);
}

void ATwiceB::serialize_a_twice_b_(std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2492, ostr);
}

void ATwiceB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_a_twice_b_(ostr);
}

}  // namespace feasst
