
#ifndef FEASST_UTILS_ARGUMENT_PARSE_H_
#define FEASST_UTILS_ARGUMENT_PARSE_H_

#include <vector>
#include <map>
#include <string>

namespace feasst {

/**
  Parse command line arguments given to main() in C++.
 */
class ArgumentParse {
 public:
  /// Construct with documentation to print when "-h" option is given.
  ArgumentParse(const char * header,
    /// Add documention and default value for options.
    /// For each option, the first value is the option itself.
    /// The second is the documentation string.
    /// And the third (optional) is the default value.
    std::vector<std::vector<std::string> > options =
    std::vector<std::vector<std::string> >());

  /// Parse the c-style input arguments.
  /// If "-h" option is given, print the documentation.
  /// Otherwise, return the arguments as a string.
  std::string parse(int argc, char ** argv);

  /// Construct directly from the c-style input arguments.
  ArgumentParse(int argc, char ** argv) { parse(argc, argv); }

  /// Return the input arguments.
  const std::vector<std::string>& args() const { return args_; }

  /// Return the input arguments as a string.
  const std::string str();

  /// Return true if option was given.
  bool option_given(const std::string& option) const;

  /// Return option.
  const std::string get(const std::string& option) const;

  /// Same as get(), but return option as a double.
  double get_double(const std::string& option) const;

  /// Same as get(), but return option as an int.
  int get_int(const std::string& option) const;

 private:
  bool is_parsed_ = false;
  std::vector<std::string> args_;
  std::vector<std::string> default_option_;
  std::vector<std::string> default_value_;
  std::string doc_;
};

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENT_PARSE_H_
