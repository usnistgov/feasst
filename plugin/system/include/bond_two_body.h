
#ifndef FEASST_CONFIGURATION_BOND_TWO_BODY_H_
#define FEASST_CONFIGURATION_BOND_TWO_BODY_H_

#include <memory>
#include "math/include/position.h"
#include "configuration/include/bond.h"

namespace feasst {

class BondTwoBody {
 public:
  BondTwoBody() {}
  virtual double energy(
      const Position& relative,
      const Bond& bond) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondTwoBody> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondTwoBody> >& deserialize_map();
  std::shared_ptr<BondTwoBody> deserialize(std::istream& istr);
  virtual ~BondTwoBody() {}

 protected:
  std::string class_name_ = "BondTwoBody";

  void serialize_bond_two_body_(std::ostream& ostr) const;
  BondTwoBody(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_TWO_BODY_H_
