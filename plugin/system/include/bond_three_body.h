
#ifndef FEASST_CONFIGURATION_BOND_THREE_BODY_H_
#define FEASST_CONFIGURATION_BOND_THREE_BODY_H_

#include "math/include/position.h"
#include "configuration/include/bond.h"

namespace feasst {

/**
  angle 0 - 1 - 2
  r01 = r0 - r1, r21 = r2 - r1
 */
class BondThreeBody {
 public:
  BondThreeBody() {}
  virtual double energy(
      const Position& relative01,
      const Position& relative21,
      const Angle& angle) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondThreeBody> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondThreeBody> >& deserialize_map();
  std::shared_ptr<BondThreeBody> deserialize(std::istream& istr);
  virtual ~BondThreeBody() {}

 protected:
  std::string class_name_ = "BondThreeBody";

  void serialize_bond_three_body_(std::ostream& ostr) const;
  BondThreeBody(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_THREE_BODY_H_
