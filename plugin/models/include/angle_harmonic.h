
#ifndef FEASST_SYSTEM_ANGLE_HARMONIC_H_
#define FEASST_SYSTEM_ANGLE_HARMONIC_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(angle) = k_energy_per_radian_sq*(angle - equilibrium_degrees)^2
  with parameters given in Angle Properties.

  The usual 1/2 factor is not included, but can be incorporated into
  the k parameter manually by the user input to the particle file.
 */
class AngleHarmonic : public BondThreeBody {
 public:
  AngleHarmonic() {}
  double energy(const double radians, const Bond& angle) const override;
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
