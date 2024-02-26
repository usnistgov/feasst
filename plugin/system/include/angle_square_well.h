
#ifndef FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_
#define FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_

#include <memory>
#include "system/include/bond_three_body.h"

namespace feasst {

/**
  U(angle) = 0 when the angle is between the minimum and maximum specified in
  AngleProperties, otherwise infinity.
  The minimum and maximum parameters are in units of degrees.
 */
class AngleSquareWell : public BondThreeBody {
 public:
  AngleSquareWell() {}
  double energy(const double radians, const Bond& angle) const override;
  std::shared_ptr<BondThreeBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AngleSquareWell(std::istream& istr);
  virtual ~AngleSquareWell() {}

 protected:
  void serialize_angle_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<AngleSquareWell> MakeAngleSquareWell() {
  return std::make_shared<AngleSquareWell>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_ANGLE_SQUARE_WELL_H_
