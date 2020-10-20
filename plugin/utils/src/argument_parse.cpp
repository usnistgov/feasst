#include <algorithm>
#include <iostream>
#include "utils/include/utils.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/argument_parse.h"

namespace feasst {

ArgumentParse::ArgumentParse(const char * header,
    std::vector<std::vector<std::string> > options) {
  doc_ = header;
  doc_.append("\n\noptions:\n");
  for (const std::vector<std::string>& opt : options) {
    doc_.append(opt[0] + ": " + opt[1]);
    if (static_cast<int>(opt.size()) == 3) {
      default_option_.push_back(opt[0]);
      doc_.append(" (default: " + opt[2] + ").");
      default_value_.push_back(opt[2]);
    }
    doc_.append("\n");
  }
}

std::string ArgumentParse::parse(int argc, char ** argv) {
  for (int index = 1; index < argc; ++index) {
    args_.push_back(std::string(argv[index]));
  }

  // find which defaults were overwritten, and use those that were not.
  std::vector<bool> overwritten_(default_option_.size(), false);
  for (int index = 1; index < argc; index += 2) {
    int findex = -1;
    if (find_in_list(std::string(argv[index]), default_option_, &findex)) {
      overwritten_[findex] = true;
    }
  }
  for (int df = 0; df < static_cast<int>(default_option_.size()); ++df) {
    if (!overwritten_[df]) {
      args_.push_back(default_option_[df]);
      args_.push_back(default_value_[df]);
    }
  }

  if (option_given("-h") || option_given("--help")) {
    std::cout << doc_ << std::endl;
    exit(0);
  }
  ASSERT(!is_parsed_, "cannot parse args twice");
  is_parsed_ = true;
  return str();
}

const std::string ArgumentParse::str() { return feasst_str(args()); }

bool ArgumentParse::option_given(const std::string& option) const {
  return std::find(args_.begin(), args_.end(), option) != args_.end();
}

const std::string ArgumentParse::get(const std::string& option) const {
  auto itr = std::find(args_.begin(), args_.end(), option);
  if (itr != args_.end() && ++itr != args_.end()) return *itr;
  FATAL("option: \"" << option << "\" not given");
}

double ArgumentParse::get_double(const std::string& option) const {
  return str_to_double(get(option));
}

int ArgumentParse::get_int(const std::string& option) const {
  return str_to_int(get(option));
}

}  // namespace feasst
