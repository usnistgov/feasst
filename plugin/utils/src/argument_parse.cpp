#include <algorithm>
#include <iostream>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/argument_parse.h"

namespace feasst {

std::string ArgumentParse::parse(int argc, char ** argv) {
  for (int index = 1; index < argc; ++index) {
    args_.push_back(std::string(argv[index]));
  }
  if (option_given("-h") || option_given("--help")) {
    std::cout << doc_ << std::endl;
    exit(0);
  }
  return str();
}

const std::string ArgumentParse::str() { return feasst_str(args()); }

bool ArgumentParse::option_given(const std::string& option) const {
  return std::find(args_.begin(), args_.end(), option) != args_.end();
}

const std::string ArgumentParse::get_(const std::string& option) const {
  auto itr = std::find(args_.begin(), args_.end(), option);
  if (itr != args_.end() && ++itr != args_.end()) return *itr;
  return std::string();
}

const std::string ArgumentParse::get(const std::string& option,
  const std::string dflt) const {
  const std::string result = get_(option);
  if (result.empty()) {
    ASSERT(!dflt.empty(), "option: " << option << " was not provided");
    return dflt;
  } else {
    return result;
  }
}

double ArgumentParse::get_double(const std::string& option) const {
  std::string arg = get_(option);
  ASSERT(!arg.empty(), "option: " << option << " was not provided");
  return str_to_double(arg);
}

double ArgumentParse::get_double(const std::string& option,
  const double dflt) const {
  std::string arg = get_(option);
  if (arg.empty()) return dflt;
  return str_to_double(arg);
}

int ArgumentParse::get_int(const std::string& option) const {
  std::string arg = get_(option);
  ASSERT(!arg.empty(), "option: " << option << " was not provided");
  return str_to_int(arg);
}

int ArgumentParse::get_int(const std::string& option,
  const int dflt) const {
  std::string arg = get_(option);
  if (arg.empty()) return dflt;
  return str_to_int(arg);
}

}
