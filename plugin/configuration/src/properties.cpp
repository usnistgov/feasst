#include "configuration/include/properties.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

void Properties::check() const {
  ASSERT(values_.size() == names_.size(),
    "size error");
  // make sure there are no spaces in property names
  for (const std::string name : names_) {
    ASSERT(num_spaces(name) == 0, "spaces are not allowed in property names");
  }
}

void Properties::add(const std::string name, const double value) {
  ASSERT(!find_in_list(name, names_),
    "property(" << name << ") already exists");
  values_.push_back(value);
  names_.push_back(name);
}

double Properties::value(const std::string name) const {
  double val;
  int index;
  ASSERT(value(name, &val, &index), "property(" << name <<") not found");
  return val;
}

bool Properties::value(const std::string name, double * val) const {
  int index;
  return value(name, val, &index);
}

bool Properties::value(const std::string name,
    double * val,
    int * index) const {
  DEBUG("finding " << name << " " << feasst_str(names_));
  bool is_found = find_in_list(name, names_, index);
  if (is_found) {
    *val = values_[*index];
  }
  return is_found;
}

void Properties::set(const std::string name, const double value) {
  int index = 0;
  const bool is_found = find_in_list(name, names_, &index);
  ASSERT(is_found, "property(" << name << ") not found");
  TRACE("setting " << name << " value " << value);
  values_[index] = value;
}

void Properties::add_or_set(const std::string name, const double value) {
  int index;
  if (find_in_list(name, names_, &index)) {
    values_[index] = value;
  } else {
    add(name, value);
  }
}

std::string Properties::str() const {
  check();
  std::stringstream ss;
  for (int index = 0; index < static_cast<int>(values_.size()); ++index) {
    ss << "{" << names_[index] << " : " << values_[index] << "}, ";
  }
  return ss.str();
}

void PropertiedEntity::set_properties(const Properties& properties,
    const std::vector<std::string>& exclude) {
  const std::vector<std::string>& name = properties_.names();
  const std::vector<double>& values = properties.values();
  ASSERT(name.size() == values.size(),
    "the same number of properties must exist to equate them");
  for (int index = 0; index < static_cast<int>(name.size()); ++index) {
    bool is_not_excl = true;
    for (const std::string& excl : exclude) {
      const std::string beginning = name[index].substr(0, excl.size());
      if (beginning == excl) {
        is_not_excl = false;
        DEBUG("excluded " << excl << " " << name[index]);
      }
    }
    if (is_not_excl) {
      set_property(index, values[index]);
    }
  }
}

bool Properties::is_equal(const Properties& properties,
    const double tolerance) const {
  return feasst::is_equal(values_, properties.values_, tolerance);
}

bool Properties::is_equal(const Properties& properties) const {
  return is_equal(properties, NEAR_ZERO);
}

void Properties::serialize(std::ostream& ostr) const {
  check();
  ostr << MAX_PRECISION;
  feasst_serialize_version(847, ostr);
  feasst_serialize(names_, ostr);
  feasst_serialize(values_, ostr);
}

Properties::Properties(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(847 == version, "version mismatch:" << version);
  feasst_deserialize(&names_, istr);
  feasst_deserialize(&values_, istr);
}

}  // namespace feasst
