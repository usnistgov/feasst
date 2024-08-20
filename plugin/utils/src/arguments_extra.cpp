#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/arguments_extra.h"

namespace feasst {

std::string str(const std::string& key, const argtype& args) {
  auto pair = args.find(key);
  if (pair != args.end()) {
    return pair->second;
  }
  return std::string();
}

std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list) {
  std::stringstream ss(line);
  std::string major;
  ss >> major;
  argtype args;
  int num_pairs = 0;
  while (!ss.eof()) {
    std::string minor, value;
    ss >> minor >> value;
    if (minor.empty()) {
      break;  // skip trailing whitespace
    }
    ASSERT(!value.empty(), "Error parsing text file on line: \"" << ss.str()
      << "\". Line syntax typically requires an odd number of space-separated "
      << "strings (e.g., Object key0 value0 ... keyN valueN. This error "
      << "typically occurs when one of a key/value pair is missing.");
    DEBUG("major " << major << " minor " << minor << " value " << value);
    if (major == "set_variable") {
      DEBUG("setting variable");
      ASSERT(variables, "set_variable found when not expected");
      (*variables)[minor] = value;
      *assign_to_list = false;
    } else if (variables && variables->count(value) > 0) {
      DEBUG("using variable");
      args[minor] = (*variables)[value];
    } else {
      DEBUG("no variable: " << value << " sz " << value.size());
      if (value.size() > 7) {
        DEBUG(value.substr(0, 7));
        if (value.substr(0, 7) == "/feasst") {
          DEBUG("replaced: " << value);
          value.replace(0, 7, install_dir());
        }
      }
      args[minor] = value;
    }
    ++num_pairs;
    ASSERT(num_pairs < 1e8, "reached maximum number of pairs");
  }
  return std::pair<std::string, argtype>(major, args);
}

argtype line_to_argtype(const std::string line) {
  argtype args;
  std::stringstream ss(line);
  std::string key, value;
  while (!ss.eof()) {
    ss >> key >> value;
    args[key] = value;
  }
  return args;
}

void replace_value(const std::string search, const std::string replace,
                   arglist * args) {
  for (std::pair<std::string, argtype>& pair1 : *args) {
    for (auto& pair2 : pair1.second) {
      if (pair2.second == search) {
        pair2.second = replace;
      }
    }
  }
}

void replace_in_value(const std::string& from, const std::string& to,
                     arglist * args) {
  for (std::pair<std::string, argtype>& pair1 : *args) {
    for (auto& pair2 : pair1.second) {
      replace(from, to, &pair2.second);
    }
  }
}

std::vector<double> parse_dimensional(const std::string& key, argtype * args,
    const int max) {
  int dim = 0;
  std::stringstream ss;
  ss << key << dim;
  std::vector<double> data;
  while (used(ss.str(), *args)) {
    data.push_back(dble(ss.str(), args));
    ASSERT(dim < max, "dim: " << dim << " is > max: " << max);
    ++dim;
    ss.str("");
    ss << key << dim;
  }
  return data;
}

}  // namespace feasst
