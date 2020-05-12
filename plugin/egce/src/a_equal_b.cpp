#include "utils/include/serialize.h"
#include "egce/include/a_equal_b.h"

namespace feasst {

AEqualB::AEqualB(const argtype& args) {
  class_name_ = "AEqualB";
  Arguments args_(args);
  extra_A_ = args_.key("extra_A").dflt("0").integer();
  num_A_ = ConstrainNumParticles({
    {"type", args_.key("particle_type_A").dflt("0").str()}});
  num_B_ = ConstrainNumParticles({
    {"type", args_.key("particle_type_B").dflt("1").str()}});
  ASSERT(num_A_.type() != num_B_.type(), "particle_type_A: " << num_A_.type()
    << " == particle_type_B_: " << num_B_.type());
}

bool AEqualB::is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const {
  const int nA = num_A_.num_particles(system, acceptance);
  const int nB = num_B_.num_particles(system, acceptance);
  bool allowed = false;
  if ( (nA == nB) || (nA == nB + extra_A_) ) {
    allowed = true;
  }
  DEBUG("nA " << nA << " nB " << nB << " allowed " << allowed);
  return allowed;
}

class MapAEqualB {
 public:
  MapAEqualB() {
    auto obj = MakeAEqualB();
    obj->deserialize_map()["AEqualB"] = obj;
  }
};

static MapAEqualB mapper_ = MapAEqualB();

std::shared_ptr<Constraint> AEqualB::create(std::istream& istr) const {
  return std::make_shared<AEqualB>(istr);
}

AEqualB::AEqualB(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "AEqualB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2492 == version, "mismatch version: " << version);
  feasst_deserialize(&extra_A_, istr);
  feasst_deserialize_fstobj(&num_A_, istr);
  feasst_deserialize_fstobj(&num_B_, istr);
}

void AEqualB::serialize_a_equal_b_(std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2492, ostr);
  feasst_serialize(extra_A_, ostr);
  feasst_serialize_fstobj(num_A_, ostr);
  feasst_serialize_fstobj(num_B_, ostr);
}

void AEqualB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_a_equal_b_(ostr);
}

}  // namespace feasst
