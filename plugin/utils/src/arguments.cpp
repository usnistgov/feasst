#include <sstream>
#include <cmath>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/io.h"

namespace feasst {

bool used(const std::string& key, const argtype& args) {
  const auto pair = args.find(key);
  if (pair != args.end()) {
    return true;
  }
  return false;
}

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

void check_all_used(const argtype& args) {
  ASSERT(args.size() == 0, "unused: " << str(args));
}

}  // namespace feasst
