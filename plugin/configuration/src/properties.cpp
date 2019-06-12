#include "configuration/include/properties.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/utils_io.h"

namespace feasst {

void Properties::check() const {
  ASSERT(property_value_.size() == property_name_.size(),
    "size error");
  // make sure there are no spaces in property names
  for (const std::string name : property_name_) {
    ASSERT(num_spaces(name) == 0, "spaces are not allowed in property names");
  }
}

void Properties::add(const std::string name, const double value) {
  ASSERT(!find_in_list(name, property_name_),
    "property(" << name << ") already exists");
  property_value_.push_back(value);
  property_name_.push_back(name);
}

double Properties::value(const std::string name) const {
  double val;
  ASSERT(value(name, &val), "property(" << name <<") not found");
  return val;
}

bool Properties::value(const std::string name, double * value) const {
  int index;
  DEBUG("finding " << name << " " << feasst_str(property_name_));
  bool found = find_in_list(name, property_name_, &index);
  if (found) {
    *value = property_value_[index];
  }
  return found;
}

void Properties::set(const std::string name, const double value) {
  int index;
  ASSERT(find_in_list(name, property_name_, &index), "property(" << name
    << ") not found");
  TRACE("setting " << name << " value " << value);
  property_value_[index] = value;
}

void Properties::add_or_set(const std::string name, const double value) {
  int index;
  if (find_in_list(name, property_name_, &index)) {
    property_value_[index] = value;
  } else {
    add(name, value);
  }
}

std::string Properties::str() const {
  check();
  std::stringstream ss;
  for (int index = 0; index < static_cast<int>(property_value_.size()); ++index) {
    ss << "{" << property_name_[index] << " : " << property_value_[index] << "} ";
  }
  return ss.str();
}

void PropertiedEntity::set_properties(const Properties& properties,
    const std::vector<std::string>& exclude) {
  const std::vector<std::string>& name = properties_.property_name();
  const std::vector<double>& values = properties.property_value();
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

void Properties::serialize(std::ostream& ostr) const {
  check();
  ostr << MAX_PRECISION;
  ostr << "1 "; // version
  feasst_serialize(property_name_, ostr);
  feasst_serialize(property_value_, ostr);
}

Properties::Properties(std::istream& istr) {
  DEBUG("deserializing properties");
  int version;
  istr >> version;
  feasst_deserialize(&property_name_, istr);
  feasst_deserialize(&property_value_, istr);
}

}  // namespace feasst
