#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"

namespace feasst {

bool used(const std::string& key, const argtype& args) {
  const auto pair = args.find(key);
  if (pair != args.end()) {
    return true;
  }
  return false;
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

float flt(const std::string& key, argtype * args,
    const float dflt) {
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

int64_t integer64(const std::string& key, argtype * args) {
  return str_to_int64(str(key, args));
}

int64_t integer64(const std::string& key, argtype * args,
    const int64_t dflt) {
  auto pair = args->find(key);
  if (pair != args->end()) {
    const std::string return_str = pair->second;
    args->erase(pair);
    return str_to_int64(return_str);
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

void feasst_check_all_used(const argtype& args) {
  if (args.size() != 0) {
    ASSERT(args.size() == 1 && args.begin()->first.empty() &&
      args.begin()->second.empty(),
      "unused argument(s): " << feasst::str(args) << ". If the arguments " <<
      "are unused, the argument parser did not expect the first keyword " <<
      "in each of the above argument pair(s). There may be a typo, a " <<
      "keyword intended for a different class, or a deprecated keyword.");
  }
}

}  // namespace feasst
