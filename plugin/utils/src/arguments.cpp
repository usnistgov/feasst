#include <sstream>
#include <cmath>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "utils/include/utils.h"
#include "utils/include/io.h"

namespace feasst {

void Arguments::init(const argtype &args) {
  args_ = args;
}

Arguments& Arguments::key(const std::string &key) {
  key_ = key;
  used_keys_.push_back(key_);
  remove_ = false;
  return *this;
}

Arguments& Arguments::dflt(const std::string &defaultVal) {
  default_value_ = defaultVal;
  return *this;
}

bool Arguments::empty() {
  ASSERT(!key_.empty(), "key must be set before");
  auto pair = args_.find(key_);
  if (pair != args_.end()) {
    return false;
  }
  return true;
}

std::string Arguments::str() {
  ASSERT(!key_.empty(), "key must be set before");
  auto pair = args_.find(key_);
  std::string str;
  if (pair != args_.end()) {
    str = pair->second;

    // remove parsed arg from args_
    if (remove_) {
      args_.erase(pair);
      remove_ = false;
    }
  } else {
    ASSERT(!default_value_.empty(), "key(" << key_ << ") is required for args "
      << "when no default value is set.");
    str = default_value_;
  }
  key_ = "";
  default_value_ = "";
  return str;
}

bool Arguments::check_all_used() {
  used_keys_.push_back("");
  for (const auto &pair : args_) {
    ASSERT(find_in_list(pair.first, used_keys_),
      "Key arg(" << pair.first <<") "
      << ") is not recognized. Check for any typos in your provided keyword "
      << "arguments and check documentation. All keywords provided in args "
      << "must be used, otherwise, a simple typo in keyword arguments would go "
      << "unnoticed. If you are a developer, check that you have processed "
      << "all possible keywords in the appropriate constructor.");
  }
  return true;
}

Arguments& Arguments::remove() {
  remove_ = true;
  return *this;
}

int str_to_int(const std::string& str) {
  std::stringstream errmsg;
  int intVal = -1;
  errmsg << str << " was " << "expected to be an integer.";
  try {
    intVal = stoi(str);
  } catch (...) {
    FATAL(errmsg.str());
  }
  const double dble = stod(str);
  ASSERT(std::abs(dble - static_cast<double>(intVal)) < 10*NEAR_ZERO,
    errmsg.str());
  return intVal;
}

bool str_to_bool(const std::string& str) {
  std::stringstream errmsg;
  errmsg << str << " was expected to be a boolean";
  if (str == "true" || str == "True") {
    return true;
  } else if (str == "false" || str == "False") {
    return false;
  } else if (str == "1") {
    return true;
  } else if (str == "0") {
    return false;
  } else {
    ASSERT(0, errmsg.str());
  }
  return -1;
}

double str_to_double(const std::string& str) {
  double double_value = -1;
  try {
    double_value = stod(str);
  } catch (...) {
    FATAL(str << " was expected to be a double precision number.");
  }
  return double_value;
}

std::string Arguments::status() const {
  std::stringstream ss;
  ss << "used(";
  for (std::string used : used_keys_) {
    ss << used << ",";
  }
  ss << "),args" << feasst_str(args_);
  return ss.str();
}

argtype Arguments::remove(const std::string first, const argtype& args) const {
  argtype new_args = args;
  auto pair = new_args.find(first);
  if (pair != new_args.end()) {
    new_args.erase(pair);
  }
  return new_args;
}

argtype Arguments::append(const std::string append,
                          const std::string first,
                          const argtype& args) const {
  argtype new_args = args;
  auto pair = new_args.find(first);
  if (pair != new_args.end()) {
    const std::string second = pair->second;
    new_args.erase(pair);
    new_args.insert({first, second + append});
  }
  return new_args;
}

//Arguments::~Arguments() {}
Arguments::~Arguments() { if (check_) check_all_used(); }

}  // namespace feasst
