
#ifndef FEASST_STEPPERS_CHIRALITY_2D_H_
#define FEASST_STEPPERS_CHIRALITY_2D_H_

#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the 2D chirality by analyzing the sign of the cross product between
  two bonds.
 */
class Chirality2D : public Analyze {
 public:
  //@{
  /** @name Arguments
    - group: group index in Configuration (default: 0).
    - bond1: index of the first bond (default: 0).
    - bond2: index of the second bond (default: 1).
    - sign_error: error if chirality of sign and != 0 (default 0).
    - Stepper arguments.
   */
  explicit Chirality2D(argtype args = argtype());
  explicit Chirality2D(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  /// Return the average number of chiral positive particles.
  const Accumulator& num_positive() const { return accumulator(); }

  // serialize
  std::string class_name() const override {
    return std::string("Chirality2D"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Chirality2D>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Chirality2D>(args); }
  Chirality2D(std::istream& istr);

  //@}
 private:
  int group_, bond1_, bond2_, sign_error_;
};

inline std::shared_ptr<Chirality2D> MakeChirality2D(
    argtype args = argtype()) {
  return std::make_shared<Chirality2D>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHIRALITY_2D_H_
