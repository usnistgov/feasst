
#include <sstream>
#include "core/include/accumulator.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

std::string Bias::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << ln_macro_prob().value(bin);
  return ss.str();
}

int Bias::bin_(
    const int macrostate_old,
    const int macrostate_new,
    const bool is_accepted) {
  int bin = -1;
  if (is_accepted) {
    bin = macrostate_new;
  } else {
    bin = macrostate_old;
  }
  DEBUG("bin " << bin);
  return bin;
}

std::map<std::string, std::shared_ptr<Bias> >& Bias::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Bias> >* ans =
     new std::map<std::string, std::shared_ptr<Bias> >();
  return *ans;
}

void Bias::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Bias> Bias::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Bias> Bias::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Bias::serialize_bias_(std::ostream& ostr) const {
  feasst_serialize_version(863, ostr);
  feasst_serialize(is_complete_, ostr);
}

Bias::Bias(std::istream& istr) {
  std::string class_name;
  istr >> class_name;
  ASSERT(!class_name.empty(), "mismatch");
  const int version = feasst_deserialize_version(istr);
  ASSERT(863 == version, "mismatch version: " << version);
  feasst_deserialize(&is_complete_, istr);
}

}  // namespace feasst
