#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "egce/include/a_equal_b.h"

namespace feasst {

AEqualB::AEqualB(argtype * args) : Constraint() {
  class_name_ = "AEqualB";
  extra_A_ = integer("extra_A", args, 0);
  num_A_ = ConstrainNumParticles({
    {"type", str("particle_type_A", args, "0")}});
  num_B_ = ConstrainNumParticles({
    {"type", str("particle_type_B", args, "1")}});
  //ASSERT(num_A_.type() != num_B_.type(), "particle_type_A: " << num_A_.type()
  //  << " == particle_type_B_: " << num_B_.type());
}
AEqualB::AEqualB(argtype args) : AEqualB(&args) {
  feasst_check_all_used(args);
}

bool AEqualB::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) {
  const int nA = num_A_.num_particles(system, acceptance);
  const int nB = num_B_.num_particles(system, acceptance);
  bool allowed = false;
  if ( (nA == nB) || (nA == nB + extra_A_) ) {
    allowed = true;
  }
  DEBUG("nA " << nA << " nB " << nB << " allowed " << allowed);
  return allowed;
}

FEASST_MAPPER(AEqualB,);

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
