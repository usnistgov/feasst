
#ifndef FEASST_CONFIGURATION_ANGLE_SQUARE_WELL_H_
#define FEASST_CONFIGURATION_ANGLE_SQUARE_WELL_H_

#include <memory>
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
      const Angle& angle) const override;
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
