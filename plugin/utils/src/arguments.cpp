#include <sstream>
#include <cmath>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
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

std::string str(const argtype& args) {
  std::stringstream ss;
  ss << "{";
  for (const auto &pair : args) {
    ss << "{" << pair.first << "," << pair.second << "},";
  }
  ss << "}";
  return ss.str();
}

}  // namespace feasst
