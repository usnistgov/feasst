#include "core/include/properties.h"
#include "core/include/debug.h"
#include "core/include/utils.h"
#include "core/include/utils_io.h"

namespace feasst {

void Properties::check_size() {
  ASSERT(property_value_.size() == property_name_.size(),
    "size error");
}

void Properties::add(const std::string name, const double value) {
  ASSERT(!find_in_list(name, property_name_), "property already exists");
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
  DEBUG("finding " << name << " " << str(property_name_));
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

}  // namespace feasst
