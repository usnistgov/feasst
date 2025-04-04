#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/bias.h"

namespace feasst {

std::string Bias::write_per_bin(const int bin) const {
  std::stringstream ss;
  ss << MAX_PRECISION << ln_prob().value(bin) << ",";
  if (bin == 0) {
    ss << "NaN";
  } else {
    ss << std::setprecision(std::cout.precision())
       << ln_prob().value(bin) - ln_prob().value(bin - 1);
  }
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

void Bias::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Bias> Bias::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Bias> Bias::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Bias> Bias::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Bias> Bias::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

void Bias::serialize_bias_(std::ostream& ostr) const {
  feasst_serialize_version(863, ostr);
  feasst_serialize(is_complete_, ostr);
  feasst_serialize(phase_, ostr);
}

Bias::Bias(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(863 == version, "mismatch version: " << version);
  feasst_deserialize(&is_complete_, istr);
  feasst_deserialize(&phase_, istr);
}

void Bias::set_cm(const int macro, const Bias& bias) {
  FATAL("not implemented");
}

const CollectionMatrix& Bias::cm() const {
  FATAL("not implemented");
}

const int Bias::visits(const int macro, const int index) const {
  FATAL("not implemented");
}

double Bias::ln_bias(const int bin_new, const int bin_old) const {
  return ln_prob().value(bin_old) - ln_prob().value(bin_new);
}

std::string Bias::write_per_bin_header(const std::string& append) const {
  std::stringstream ss;
  ss << "ln_prob" << append << ",delta_ln_prob" << append;
  return ss.str();
}

}  // namespace feasst
