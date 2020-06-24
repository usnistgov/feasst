
#ifndef FEASST_MATH_SOLVER_H_
#define FEASST_MATH_SOLVER_H_

#include <map>
#include <memory>
#include <string>
#include <sstream>
#include "utils/include/arguments.h"

namespace feasst {

class Formula;

/**
  Solver is used to find roots.
 */
class Solver {
 public:
  /**
    args:
    - tolerance: solve within this tolerance.
    - lower: optional lower bound.
    - upper: optional upper bound.
    - guess: initial guess for the root.
   */
  Solver(const argtype& args = argtype());

  /// Return the tolerance
  double tolerance() const { return tolerance_; }

  /// Return true if lower bound is used
  bool is_lower() const { return is_lower_; }

  /// Return lower bound
  double lower() const { return lower_; }

  /// Set lower bound
  void set_lower(const double lower) { is_lower_ = true; lower_ = lower; }

  /// Return true if upper bound is used
  bool is_upper() const { return is_upper_; }

  /// Return upper bound
  double upper() const { return upper_; }

  /// Set upper bound
  void set_upper(const double upper) { is_upper_ = true; upper_ = upper; }

  /// Return the guess
  double guess() const { return guess_; }

  /// Set the guess
  void set_guess(const double guess) { is_guess_ = true; guess_ = guess; }

  /// Return a root of the Formula.
  virtual double root(Formula * formula) = 0;

  // serialization
  virtual std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Solver> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Solver> >& deserialize_map();
  std::shared_ptr<Solver> deserialize(std::istream& istr);
  explicit Solver(std::istream& istr);
  virtual ~Solver() {}

 protected:
  std::string class_name_ = "Solver";
  void serialize_solver_(std::ostream& ostr) const;
  Arguments args_;

 private:
  double tolerance_;
  double lower_ = 0;
  double upper_ = 0;
  double guess_ = 0;
  bool is_lower_ = false;
  bool is_upper_ = false;
  bool is_guess_ = false;
};

}  // namespace feasst

#endif  // FEASST_MATH_SOLVER_H_
