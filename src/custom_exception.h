/**
 * Authors: Nathan Mahynksi wrote this class
 *          Harold Hatch editted
 */
#ifndef CUSTOM_EXCEPTION_H_
#define CUSTOM_EXCEPTION_H_

#include <exception>
#include <string>
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;

namespace feasst {

/**
 * Customize the way that exceptions are handled and thrown.
 */
class CustomException : public std::exception {
 public:
  /// Instantiate an exeption with a user-defined error message.
  CustomException(std::string m = "custom exception occurred");

  /// Instantiate an exeption with a user-defined error message.
  CustomException(std::stringstream& m);

  ~CustomException() throw() { cout << msg_; }

  /// Add additional information to the message, such as the thread number.
  void catMessage();

  /// Return the user's message
  const char* what() const throw() { return msg_.c_str(); }

 protected:
  std::string msg_;
};

}  // namespace feasst

#endif  // CUSTOM_EXCEPTION_H_
