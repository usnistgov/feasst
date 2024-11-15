
#ifndef FEASST_CHAIN_ANALYZE_BONDS_H_
#define FEASST_CHAIN_ANALYZE_BONDS_H_

#include <vector>
#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "system/include/rigid_angle.h"
#include "system/include/rigid_dihedral.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the distribution and moments of each type of bond and angle.
 */
class AnalyzeBonds : public Analyze {
 public:
  //@{
  /** @name Arguments
    - bond_bin_width: bin width for bond histogram (default: 1).
    - bond_bin_center: center of first bin in bond histogram (default: 0).
    - angle_bin_width: bin width for angle histogram in units of radians.
      (default: 0.01).
    - angle_bin_center: center of first bin in angle histogram (default: 0).
    - dihedral_bin_width: bin width for dihedral histogram in units of radians.
      (default: 0.01).
    - dihedral_bin_center: center of first bin in dihedral histogram (default: 0).
    - Stepper arguments.
   */
  explicit AnalyzeBonds(argtype args = argtype());
  explicit AnalyzeBonds(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(const MonteCarlo& mc) override;
  void initialize(MonteCarlo * mc) override;

  /// Return the average bond length.
  const Accumulator& bond(const int type) const { return bond_[type]; }

  /// Return the average angle in radians.
  const Accumulator& angle(const int type) const { return angle_[type]; }

  /// Return the average dihedral in radians.
  const Accumulator& dihedral(const int type) const { return dihedral_[type]; }

  /// Return the histogram of bond length.
  const Histogram& bond_hist(const int type) const { return bond_hist_[type]; }

  /// Return the histogram of angle in radians.
  const Histogram& angle_hist(const int type) const {
    return angle_hist_[type]; }

  /// Return the histogram of dihedral in radians.
  const Histogram& dihedral_hist(const int type) const {
    return dihedral_hist_[type]; }

  std::string write(const MonteCarlo& mc) override;

  std::string class_name() const override { return std::string("AnalyzeBonds"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeBonds>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<AnalyzeBonds>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeBonds(std::istream& istr);

  //@}
 private:
  std::vector<Accumulator> bond_, angle_, dihedral_;
  std::vector<Histogram> bond_hist_, angle_hist_, dihedral_hist_;

  // temporary and not serialized
  RigidAngle angle_calc_;
  RigidDihedral dihedral_calc_;
};

inline std::shared_ptr<AnalyzeBonds> MakeAnalyzeBonds(
    argtype args = argtype()) {
  return std::make_shared<AnalyzeBonds>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_ANALYZE_BONDS_H_
