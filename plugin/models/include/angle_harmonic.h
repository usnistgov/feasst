
#ifndef FEASST_SYSTEM_ANGLE_HARMONIC_H_
#define FEASST_SYSTEM_ANGLE_HARMONIC_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  This class uses the following Angle Properties
  - k_energy_per_radian_sq
  - equilibrium_degrees
  - num_jacobian_gaussian (default: 1)
  which are described as follows.

  U(angle) = k_energy_per_radian_sq*(angle - equilibrium_degrees)^2
  with parameters given in Angle Properties.

  The usual 1/2 factor is not included, but can be incorporated into
  the k parameter manually by the user input to the particle file.

  The Jacobian-Gaussian algorithm described in:

  https://doi.org/10.1021/acs.jctc.7b00173

  is implemented here for the "2-branch" case, while "3-branch" is not currently
  implemented.
  Thank you to Dr Daniel Siderius and Professor Bin Chen for providing example
  implementation code and discussion.

  The num_jacobian_gaussian parameter is only read by the one angle that has
  both sites that are placed by the branch in the vertices.
  To use the default BondThreeBody::random_branch,
  set num_jacobian_gaussian to 0.
 */
class AngleHarmonic : public BondThreeBody {
 public:
  AngleHarmonic() {}
  double energy(const double radians, const Bond& angle) const override;
// HWH tests on propane revealed this to slow down the simulation.
//  double random_angle_radians(const Angle& angle, const double beta,
//    const int dimension, Random * random) const override;
  void random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    const bool is_position_held,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random,
    double * ln_met,
    const Position * const a1,
    const Position * const a2,
    const Position * const m1,
    const Position * const m2) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AngleHarmonic(std::istream& istr);
  virtual ~AngleHarmonic() {}

 protected:
  void serialize_angle_harmonic_(std::ostream& ostr) const;
};

inline std::shared_ptr<AngleHarmonic> MakeAngleHarmonic() {
  return std::make_shared<AngleHarmonic>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ANGLE_HARMONIC_H_
