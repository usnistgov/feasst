/**
 * \file
 *
 * \brief custom exception class
 *
 * Authors: Nathan Mahynksi started this class
 *          Harold Hatch editted
 */

#ifndef EXCEPTION_H_
#define EXCEPTION_H_

#include <exception>
#include <string>
#include <iostream>
#include <sstream>

using std::cout;
using std::endl;

class customException : public std::exception {
 public:

  /// Instantiate an exeption with a user-defined error message
  customException(std::string m = "custom exception occurred");
  customException(std::stringstream& m);
  ~customException();

  /// add additional information to the message
  void catMessage();

  //!< Return the user's message
  const char* what() const throw() {return msg_.c_str();}

 protected:
  std::string msg_;
};

#endif  // EXCEPTION_H_
