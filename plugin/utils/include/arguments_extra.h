
#ifndef FEASST_UTILS_ARGUMENTS_EXTRA_H_
#define FEASST_UTILS_ARGUMENTS_EXTRA_H_

#include <utility>
#include <string>
#include <vector>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Read an argument, but do not remove it.
  WARNING: This should typically only be used for error checking.
 */
std::string str(const std::string& key, const argtype& args);

/// Parse a text interface line.
/// These typically begin with an object Name, then space-separated arguments.
/// Each space seprated argument begins with the argument name, and equal sign
/// and then the argument values.
/// Backwards compatible with the old interface which uses a space instead of
/// an equal sign.
/// First, look for set_variable, to generate a list of variables to use for
/// name substitution.
/// Finally, replace any value beginnig with /feasst with the install_dir().
std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list);

/// convert a space-delimited string into argtype
argtype line_to_argtype(const std::string line);

/// Parse MonteCarlo arguments from a stream
std::vector<arglist> parse_mcs(std::istream& is, argtype variables = argtype());

/// Find all values equal to "search" in args and replace with "replace"
void replace_value(const std::string search, const std::string replace,
                   arglist * args);

void replace_in_value(const std::string& from, const std::string& to,
                      arglist * args);

/// Read data from arguments beginning with key and counting from 0 up.
/// For example, "{{"x0", "1"}, {"x1", "2"}}" will return {1, 2} vector.
std::vector<double> parse_dimensional(const std::string& key, argtype * args,
  const int max);

std::vector<std::pair<std::string, std::vector<std::string> > > parse_for(const argtype& args);

arglist expand_for(const arglist& list, int * first_end_for, int *last_for);

arglist parse_if(const arglist& list, int * first_end_if, int *last_if, int *last_else);

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENTS_EXTRA_H_
