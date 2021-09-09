
#ifndef FEASST_SYSTEM_BOND_SQUARE_WELL_H_
#define FEASST_SYSTEM_BOND_SQUARE_WELL_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  U(r) = 0 when distance is between the minimum and maximum specified in
  BondProperties, otherwise infinity.
 */
class BondSquareWell : public BondTwoBody {
 public:
  BondSquareWell() {}
  double energy(const double distance, const Bond& bond) const override;
  double random_distance(const Bond& bond, const double beta, const int dimen,
    Random * random) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit BondSquareWell(std::istream& istr);
  virtual ~BondSquareWell() {}

 protected:
  void serialize_bond_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<BondSquareWell> MakeBondSquareWell() {
  return std::make_shared<BondSquareWell>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_BOND_SQUARE_WELL_H_
