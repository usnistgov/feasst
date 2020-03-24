
#ifndef FEASST_CONFIGURATION_BOND_SQUARE_WELL_H_
#define FEASST_CONFIGURATION_BOND_SQUARE_WELL_H_

#include <memory>
#include "system/include/bond_two_body.h"

namespace feasst {

/**
  U(r) = 0 when |l-l0| < delta/2, otherwise infinity.
 */
class BondSquareWell : public BondTwoBody {
 public:
  explicit BondSquareWell(const argtype& args = argtype()) {}
  double energy(
      const Position& relative,
      const Bond& bond) const override;
  std::shared_ptr<BondTwoBody> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit BondSquareWell(std::istream& istr);
  virtual ~BondSquareWell() {}

 protected:
  void serialize_bond_square_well_(std::ostream& ostr) const;
};

inline std::shared_ptr<BondSquareWell> MakeBondSquareWell(
    const argtype &args = argtype()) {
  return std::make_shared<BondSquareWell>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_SQUARE_WELL_H_
