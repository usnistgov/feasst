#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "egce/include/a_half_b.h"

namespace feasst {

AHalfB::AHalfB(argtype * args) : Constraint() {
  class_name_ = "AHalfB";
  extra_ = integer("extra", args, 0);
  num_A_ = ConstrainNumParticles({
    {"type", str("particle_type_A", args, "0")}});
  num_B_ = ConstrainNumParticles({
    {"type", str("particle_type_B", args, "1")}});
  //ASSERT(num_A_.type() != num_B_.type(), "particle_type_A: " << num_A_.type()
  //  << " == particle_type_B_: " << num_B_.type());
}
AHalfB::AHalfB(argtype args) : AHalfB(&args) {
  feasst_check_all_used(args);
}

bool AHalfB::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) {
  const int nA = num_A_.num_particles(system, acceptance);
  const int nB = num_B_.num_particles(system, acceptance);
  bool allowed = false;
  if ( std::abs(2*nA - nB) <= extra_) {
    allowed = true;
  }
  DEBUG("nA " << nA << " nB " << nB << " extra " << extra_
    << " allowed " << allowed);
  return allowed;
}

FEASST_MAPPER(AHalfB,);

AHalfB::AHalfB(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "AHalfB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2492 == version, "mismatch version: " << version);
  feasst_deserialize(&extra_, istr);
  feasst_deserialize_fstobj(&num_A_, istr);
  feasst_deserialize_fstobj(&num_B_, istr);
}

void AHalfB::serialize_a_half_b_(std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2492, ostr);
  feasst_serialize(extra_, ostr);
  feasst_serialize_fstobj(num_A_, ostr);
  feasst_serialize_fstobj(num_B_, ostr);
}

void AHalfB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_a_half_b_(ostr);
}

}  // namespace feasst
