
#ifndef FEASST_UTILS_ARGUMENTS_H_
#define FEASST_UTILS_ARGUMENTS_H_

#include <vector>
#include <map>
#include <string>

namespace feasst {

/// Use a map of string pairs as a dictionary for arguments.
typedef std::map<std::string, std::string> argtype;
typedef std::map<std::string, argtype> arglist;

/**
 * The Arguments class takes in arguments as a map of strings
 * (e.g., key,value pairs in a dictionary).
 * It then uses chainsetters to access the key values and finds keys that
 * were used for error checking.
 */
class Arguments {
 public:
  Arguments() {}

  /**
    Initialize the arguments using a brace enclosed initializer list.
   */
  void init(const argtype &args);

  /// Construct via brace enclosed initializer list.
  explicit Arguments(const argtype &args) { init(args); }

  /**
   * Set the argument key (the first in the map/dictionary pair).
   * Store the key for error processing.
   * Return self for chainsetting.
   */
  Arguments& key(const std::string &key);

  /// Set the default value if key is not present in args.
  Arguments& dflt(const std::string &defaultVal);

  /// Return true if key is not present in args. Otherwise, false.
  bool empty();

  /// Return true if the key is used (e.g., the inverse of empty()).
  bool used() { return !empty(); }

  /// Return the value of the processed keyword.
  /// Reset key and dflt
  std::string str();

  /// Return the conversion of a str of the processed keyword to double
  /// a precision floating point number.
  double dble();

  /// Return the conversion of a str of the processed keyword to int.
  int integer();

  /// Return the conversion of a str of the processed keyword to boolean.
  /// Accept the strings "True", "False", "1", "0", "true" or "false".
  bool boolean();

  /// Upon destruction, check that all provided args were processed.
  /// Automatically include empty string key as processed.
  bool check_all_used();

  /// Chainset argparse to remove the next arg returned by str from args
  Arguments& remove();

  /// Return the size of the args (e.g., the number)
  int size() const { return static_cast<int>(args_.size()); }

  argtype args() const { return args_; }  //!< Return args

  /// Print the status of the arguments to human readable string.
  std::string status() const;

  /// Don't check for unused arguments upon destruction.
  void dont_check() { check_ = false; }

  // The following functions ease manipulation of argtypes

  /// Return argtype after removing any pair with the given first index.
  argtype remove(const std::string first, const argtype& args) const;

  /// Return argtype after appending onto pair with given first index.
  argtype append(const std::string append,
                 const std::string first,
                 const argtype& args) const;

  ~Arguments() { if (check_) check_all_used(); }

 private:
  argtype args_;
  std::string key_;
  std::string default_value_;
  std::vector<std::string> used_keys_;
  bool remove_ = false;
  bool check_ = true;
};

}  // namespace feasst

#endif  // FEASST_UTILS_ARGUMENTS_H_
