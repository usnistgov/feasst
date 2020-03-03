
#ifndef FEASST_CONFIGURATION_ANGLE_SQUARE_WELL_H_
#define FEASST_CONFIGURATION_ANGLE_SQUARE_WELL_H_

#include <math.h>
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(theta) = 0 when |theta-theta0| < delta/2, otherwise infinity.
  where theta is in degrees
 */
class AngleSquareWell : public BondThreeBody {
 public:
  explicit AngleSquareWell(const argtype& args = argtype()) {}
  double energy(
      const Position& relative01,
      const Position& relative21,
      const Angle& angle) const override {
    const double theta0 = angle.property("theta0");
    const double delta = angle.property("delta");
    const double theta = radians_to_degrees(acos(relative01.cosine(relative21)));
    TRACE("theta " << theta);
    if (std::abs(theta - theta0) > 0.5*delta) {
      return NEAR_INFINITY;
    }
    return 0.;
  }
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AngleSquareWell(std::istream& istr);
  virtual ~AngleSquareWell() {}

 protected:
  void serialize_angle_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<AngleSquareWell> MakeAngleSquareWell(
    const argtype &args = argtype()) {
  return std::make_shared<AngleSquareWell>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_ANGLE_SQUARE_WELL_H_
