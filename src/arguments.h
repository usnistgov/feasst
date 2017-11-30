/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ARGUMENTS_H_
#define ARGUMENTS_H_

#include <vector>
#include <iostream>
#include "functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * The Arguments class takes in arguments as a map of strings
 * (e.g., key,value pairs in a dictionary).
 * It then uses chainsetters to access the key values and finds keys that
 * were used for error checking.
 */
class Arguments {
 public:
  typedef std::map<std::string, std::string> argtype;
  void initArgs(const std::string className, const argtype &args) {
    className_ = className;
    args_ = args;
  }

  /**
   * Set the argument key (the first in the map/dictionary pair).
   * Store the key for error processing.
   * Return self for chainsetting.
   */
  Arguments& key(const std::string &key) {
    key_ = key;
    usedKeys_.push_back(key_);
    return *this;
  }

  /// Set the default value if key is not present in args.
  Arguments& dflt(const std::string &defaultVal) {
    defaultVal_ = defaultVal;
    return *this;
  }

  /// Return the value of the processed keyword.
  /// Reset key and dflt
  std::string str() {
    ASSERT(!key_.empty(), "key must be set before str");
    auto pair = args_.find(key_);
    std::string str;
    if (pair != args_.end()) {
      str = pair->second;
    } else {
      ASSERT(!defaultVal_.empty(), "key(" << key_ << ") is required for args "
        << "when no default value is set.");
      str = defaultVal_;
    }
    key_ = "";
    defaultVal_ = "";
    return str;
  }

  /// Upon destruction, check that all provided args were processed.
  /// Automatically include empty string key as processed.
  bool checkAllArgsUsed(){
    usedKeys_.push_back("");
    for (const auto &pair : args_) {
      ASSERT(findInList(pair.first, usedKeys_),
        "Key arg(" << pair.first <<") given to class(" << className_
        << ") is not recognized. All keywords provided in args must be "
        << "used, otherwise, a simple typo in keyword arguments would go "
        << "unnoticed. If you are a developer, check that you have processed "
        << "all possible keywords in the appropriate constructor. If you are "
        << "a user, check for any typos in your provided keyword arguments and "
        << "check documentation for the class(" << className_ << ")");
    }
    return true;
  }

 private:
  std::string className_;
  argtype args_;
  std::string key_;
  std::string defaultVal_;
  std::vector<std::string> usedKeys_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ARGUMENTS_H_
