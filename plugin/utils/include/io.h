
#ifndef FEASST_UTILS_IO_H_
#define FEASST_UTILS_IO_H_

#include <string>
#include <sstream>
#include <vector>
#include <deque>
#include <map>
#include <cstdint>
#include "utils/include/max_precision.h"
#include "utils/include/definitions.h"

namespace feasst {

/// Return string representation of vector
template<class T>
std::string feasst_str(const std::vector<T> &vec,
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
    ss << vec[i];
    if (i != static_cast<int>(vec.size()-1)) ss << ",";
  }
  return ss.str();
}
template<class T>
std::string feasst_str(const std::vector<std::vector<T> > &vec,
  const bool max_precision = false) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << feasst_str(vec[i], max_precision) << std::endl;
  }
  return ss.str();
}
template<class T>
std::string feasst_str(const std::deque<T> &deq,
  const bool max_precision = false) {
  std::stringstream ss;
  if (max_precision) ss << MAX_PRECISION;
  for (unsigned int i = 0; i < deq.size(); ++i) {
    ss << deq[i] << ",";
  }
  return ss.str();
}

/// Print a std::map of a pair of strings in human readable format.
std::string str(const std::map<std::string, std::string>& mp);

/// Split a string into a vector of strings via given delimitor.
std::vector<std::string> split(const std::string& str_to_split,
  const char delimitor = ' ');

/// Return phrase with all characters up to the last specialchr removed.
std::string trim(const char* specialchr, const char* phrase,
  /** If 1, remove characters from the left. Otherwise, remove characters
   *  from the right.*/
  int from_left = 1);

/// Same as above except for string phrase.
std::string trim(const char* specialchr, std::string phrase, int from_left = 1);

/// Convert to string with maximum precision
template <typename T>
std::string str(const T a_value) {
  std::ostringstream out;
  out << MAX_PRECISION << a_value;
  return out.str();
}

/// Return the number of spaces in string
int num_spaces(const std::string str);

/// Return true if substring is found inside of string.
bool is_found_in(const std::string& str, const std::string& substr);

/// Replace substring "from" in "str" to "to". Return true if replaced.
/// https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
bool replace(const std::string& from, const std::string& to, std::string * str);

/// Convert a string to an integer.
int str_to_int(const std::string& str,
  /// Return a FATAL error if impossible. Otherwise, return -1.
  const bool fatal = true);

/// Convert a string to a 64-bit integer.
int64_t str_to_int64(const std::string& str,
  /// Return a FATAL error if impossible. Otherwise, return -1.
  const bool fatal = true);

/// Convert a string to an double.
double str_to_double(const std::string& str);

/// Convert a string to a boolean.
bool str_to_bool(const std::string& str);

/// Print an integer with the maximum number of digits as given largest integer.
std::string sized_int_to_str(const int num, const int largest_num);

}  // namespace feasst

#endif  // FEASST_UTILS_IO_H_
