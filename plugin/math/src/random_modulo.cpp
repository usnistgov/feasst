#include <limits>
#include "math/include/random_modulo.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

double RandomModulo::gen_uniform_() {
  TRACE(std::numeric_limits<int>::max());
  TRACE("M " << M_);
  TRACE("a " << a_)
  TRACE("X " << X_)
  //X_ = (a_*X_) % M_;
  X_ = a_*X_ - M_*int(a_*X_/M_);
  TRACE("X " << X_)
  return static_cast<double>(X_)/static_cast<double>(M_);
}

void RandomModulo::reseed_(const int seed) {
  TRACE("seed " << seed << " address " << this);
  X_ = seed;
}

class MapRandomModulo {
 public:
  MapRandomModulo() {
    RandomModulo().deserialize_map()["RandomModulo"] = MakeRandomModulo();
  }
};

static MapRandomModulo mapper_ = MapRandomModulo();

RandomModulo::RandomModulo(const argtype& args) : Random(args) {
  class_name_ = "RandomModulo";
  M_ = std::pow(2, 31) - 1;
  a_ = std::pow(7, 5);
  parse_seed_(args);
}

void RandomModulo::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_random_(ostr);
  feasst_serialize_version(5935, ostr);
  feasst_serialize(M_, ostr);
  feasst_serialize(a_, ostr);
  feasst_serialize(X_, ostr);
  feasst_serialize_endcap("RandomModulo", ostr);
}

RandomModulo::RandomModulo(std::istream& istr)
  : Random(istr) {
  ASSERT(class_name_ == "RandomModulo", "name: " << class_name_);
  const int verison = feasst_deserialize_version(istr);
  ASSERT(verison == 5935, "version");
  feasst_deserialize(&M_, istr);
  feasst_deserialize(&a_, istr);
  feasst_deserialize(&X_, istr);
  feasst_deserialize_endcap("RandomModulo", istr);
}

}  // namespace feasst
