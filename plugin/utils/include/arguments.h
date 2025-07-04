
#ifndef FEASST_UTILS_ARGUMENTS_H_
#define FEASST_UTILS_ARGUMENTS_H_

#include <utility>
#include <vector>
#include <map>
#include <string>
#include <cstdint>

namespace feasst {

/// Use a map of string pairs as a dictionary for arguments.
typedef std::map<std::string, std::string> argtype;
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

/// Return true if key is used in args.
bool used(const std::string& key, const argtype& args);

/// Read an argument and remove it
std::string str(const std::string& key, argtype * args);

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

/// Read an argument and remove it from args, then return as 64-bit integer.
int64_t integer64(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
int64_t integer64(const std::string& key, argtype * args,
  const int64_t dflt);

/// Read an argument and remove it from args, then return as a boolean.
bool boolean(const std::string& key, argtype * args);

/// Same as above, but with a default value should key not be in args.
bool boolean(const std::string& key, argtype * args,
  const bool dflt);

/// Check that all arguments are used.
void feasst_check_all_used(const argtype& args);

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENTS_H_
