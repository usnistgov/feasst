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

/// Use a map of string pairs to function as a dictionary for arguments.
typedef std::map<std::string, std::string> argtype;

/**
 * The Arguments class takes in arguments as a map of strings
 * (e.g., key,value pairs in a dictionary).
 * It then uses chainsetters to access the key values and finds keys that
 * were used for error checking.
 */
class Arguments {
 public:
  void initArgs(const std::string className, const argtype &args);

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

  /// Return the value of the processed keyword.
  /// Reset key and dflt
  std::string str();

  /// Upon destruction, check that all provided args were processed.
  /// Automatically include empty string key as processed.
  bool checkAllArgsUsed();

  /// Chainset argparse to remove the next arg returned by str from args
  Arguments& rm();

  /// Return the size of the args (e.g., the number)
  int size() const { return static_cast<int>(args_.size()); }

  argtype args() const { return args_; }  //!< Return args

 private:
  std::string className_;
  argtype args_;
  std::string key_;
  std::string defaultVal_;
  std::vector<std::string> usedKeys_;
  bool removeArg_ = false;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ARGUMENTS_H_
