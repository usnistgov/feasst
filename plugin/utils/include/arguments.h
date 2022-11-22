
#ifndef FEASST_UTILS_ARGUMENTS_H_
#define FEASST_UTILS_ARGUMENTS_H_

#include <vector>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/debug.h"

namespace feasst {

/// Use a map of string pairs as a dictionary for arguments.
typedef std::map<std::string, std::string> argtype;
//typedef std::map<std::string, argtype> arglist;
typedef std::vector<std::pair<std::string, argtype> > arglist;

/**
  Classes should take argtype as input for (optional) arguments that may or may
  not have default values:

  MakeTestArgs({{"arg", "value"}});

  All arguments are pairs of strings that may be converted to int, double or
  bool via utils/include/io.h.

  Constructurs should take argtype as object and pointer.
  The point implementation removes arguments as they are processed, while the
  argtype constructor passes the arguments to the other constructor, and then
  checks that all argtype were used.
  Thus, the pointer implementation is used for base classes and containers,
  while the object implementation is for the user.

  Please see utils/test/arguments.cpp TestArgs and TEST(Arguments, args)
  for a working example of how to implement arguments in a class.
 */

/// Return true if key is used in args (templated for arglist or argtype).
template <typename T>
bool used(const std::string& key, const T& args) {
  const auto pair = args.find(key);
  if (pair != args.end()) {
    return true;
  }
  return false;
}

/**
  Read an argument, but do not remove it.
  WARNING: This should typically only be used for error checking.
 */
std::string str(const std::string& key, const argtype& args);

/// Read an argument and remove it
std::string str(const std::string& key, argtype * args);

//// Depreciate
//argtype get(const std::string& key, arglist * args);

/// Return a human-readable string representing argtype (use io template)
//std::string str(const argtype& args);

/// Return a human-readable string representing arglist
std::string str(const arglist& args);

/// Same as above, but with a default value should key not be in args.
std::string str(const std::string& key, argtype * args,
  const std::string dflt);

/// Read an argument and remove it from args, then return as double.
double dble(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
double dble(const std::string& key, argtype * args,  const double dflt);

/// Same as above, but with a default value should key not be in args.
float flt(const std::string& key, argtype * args, const float dflt);

/// Read an argument and remove it from args, then return as double.
int integer(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
int integer(const std::string& key, argtype * args,
  const int dflt);

/// Read an argument and remove it from args, then return as a boolean.
bool boolean(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
bool boolean(const std::string& key, argtype * args,
  const bool dflt);

/// Append to given key
void append(const std::string& key, argtype * args, const std::string& append);

/// Check that all arguments are used.
# define FEASST_CHECK_ALL_USED(args) \
{ \
  if (args.size() != 0) { \
    ASSERT(args.size() == 1 && args.begin()->first.empty() && \
      args.begin()->second.empty(), \
      "unused argument(s): " << feasst::str(args) << ". If the arguments are unused " << \
      "then that means the objects did not expect the first keyword " << \
      "supplied in the above argument pair(s). Thus, there was likely a " << \
      "typo or the keyword is intended for a different class."); \
  } \
}

/// If args contains derived class of T, return factory pointer and remove from
/// args.
template <class T>
std::shared_ptr<T> parse(T * obj, arglist * args) {
  std::shared_ptr<T> new_obj;
  const auto& map = obj->deserialize_map();
  //for (auto& arg : *args) {
  //for (const auto& mp : map) {
  //for (int iarg = 0; iarg < static_cast<int>(args->size()); ++iarg) {
//    int find = -1;
//  int iarg = 0;
  //INFO("parsing " << args->begin()->first);
  //for (const auto& mp : map) INFO(mp.first);
  if (map.count(args->begin()->first) > 0) {
  //if (map.count((*args)[iarg].first) > 0) {
  //if (find_in_map((*args)[iarg].first, map, &find)) {
    new_obj = obj->factory(args->begin()->first, &args->begin()->second);
    //INFO(new_obj->class_name());
    FEASST_CHECK_ALL_USED(args->begin()->second);
    //new_obj = obj->factory((*args)[iarg].first, &(*args)[iarg].second);
    args->erase(args->begin());
    //args->erase(args->begin() + iarg);
  //auto pair = args->find(mp.first);//map.find(arg.first);
  //if (pair != args->end()) {
  //  new_obj = obj->factory(pair->first, &pair->second);
  //  args->erase(pair);

    return new_obj;
  }
  //}
  return new_obj;
}

/// If an argument is not used, add it.
void add_if_not_used(const std::string& key, argtype * args,
  const std::string& value);

/// convert a space-delimited string into argtype
argtype line_to_argtype(const std::string line);

/// Find all values equal to "search" in args and replace with "replace"
void replace_value(const std::string search, const std::string replace,
                   arglist * args);

/// Find all values that contain "search" in args and replace with "replace"
void replace_in_value(const std::string& from, const std::string& to,
                      arglist * args);

/// Read data from arguments beginning with key and counting from 0 up.
/// For example, "{{"x0", "1"}, {"x1", "2"}}" will return {1, 2} vector.
std::vector<double> parse_dimensional(const std::string& key, argtype * args, const int max);

/// Parse a text interface line.
/// These typically begin with an object Name, then space-separated arguments.
/// First, look for set_variable, to generate a list of variables to use for
/// name substitution.
/// Finally, replace any value beginnig with /feasst with the install_dir().
std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list);

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENTS_H_
