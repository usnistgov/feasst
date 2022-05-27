#include <sstream>
#include <cmath>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/io.h"

namespace feasst {

std::string str(const std::string& key, const argtype& args) {
  auto pair = args.find(key);
  if (pair != args.end()) {
    return pair->second;
  }
  return std::string();
}

std::string str(const std::string& key, argtype * args) {
  auto pair = args->find(key);
  ASSERT(pair != args->end(), "key(" << key << ") is required for args but " <<
    "not found in " << str(*args));
  const std::string second = pair->second;
  args->erase(pair);
  return second;
}

//argtype get(const std::string& key, arglist * args) {
//  auto pair = args->find(key);
//  ASSERT(pair != args->end(), "key(" << key << ") is required for args but " <<
//    "not found");// << str(*args));
//  const argtype second = pair->second;
//  args->erase(pair);
//  return second;
//}

std::string str(const std::string& key, argtype * args,
    const std::string dflt) {
  std::string return_str;
  auto pair = args->find(key);
  if (pair != args->end()) {
    return_str = pair->second;
    args->erase(pair);
  } else {
    return_str = dflt;
  }
  return return_str;
}

double dble(const std::string& key, argtype * args) {
  return str_to_double(str(key, args));
}

double dble(const std::string& key, argtype * args,
    const double dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_double(return_str);
  } else {
    return dflt;
  }
}

int integer(const std::string& key, argtype * args) {
  return str_to_int(str(key, args));
}

int integer(const std::string& key, argtype * args,
    const int dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_int(return_str);
  } else {
    return dflt;
  }
}

bool boolean(const std::string& key, argtype * args) {
  return str_to_bool(str(key, args));
}

bool boolean(const std::string& key, argtype * args,
    const bool dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_bool(return_str);
  } else {
    return dflt;
  }
}

void append(const std::string& key, argtype * args, const std::string& append) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string second = pair->second;
    args->erase(pair);
    args->insert({key, second + append});
  }
}

std::string str(const arglist& args) {
  std::stringstream out;
  out << "{{";
  for (const auto& arg : args) {
    out << "{\"" << arg.first << "\",";
    out << str(arg.second) << "},";
  }
  out << "}}";
  return out.str();
}

void add_if_not_used(const std::string& key, argtype * args,
  const std::string& value) {
  if (!used(key, *args)) {
    args->insert({key, value});
  }
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

std::vector<double> parse_dimensional(const std::string& key, argtype * args, const int max) {
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

std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list) {
  std::stringstream ss(line);
  std::string major;
  ss >> major;
  argtype args;
  while(!ss.eof()) {
    std::string minor, value;
    ss >> minor >> value;
    ASSERT(!value.empty(), "Error parsing text file on line: \"" << ss.str()
      << "\". Line syntax typically requires an odd number of space-separated "
      << "strings (e.g., Object key0 value0 ... keyN valueN."
      << " This error typically occurs when one of a key/value pair is missing.");
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
  }
  return std::pair<std::string, argtype>(major, args);
}

}  // namespace feasst
