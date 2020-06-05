
#ifndef FEASST_CONFIGURATION_BOND_VISITOR_H_
#define FEASST_CONFIGURATION_BOND_VISITOR_H_

#include <memory>
#include <string>
#include <map>
#include "configuration/include/bond.h"
#include "system/include/bond_two_body.h"
#include "system/include/bond_three_body.h"

namespace feasst {

class Configuration;
class Select;

class BondVisitor {
 public:
  explicit BondVisitor(const argtype& args = argtype()) {}
  void compute(
      const BondTwoBody& model,
      const Configuration& config,
      const int group_index = 0);
  void compute(
      const BondTwoBody& model,
      const Select& selection,
      const Configuration& config);
  void compute(
      const BondThreeBody& model,
      const Configuration& config,
      const int group_index = 0);
  void compute(
      const BondThreeBody& model,
      const Select& selection,
      const Configuration& config);
  void set_energy(const double energy) { energy_ = energy; }
  double energy() const { return energy_; }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<BondVisitor> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<BondVisitor> >& deserialize_map();
  std::shared_ptr<BondVisitor> deserialize(std::istream& istr);
  explicit BondVisitor(std::istream& istr);
  virtual ~BondVisitor() {}

 protected:
  std::string class_name_ = "BondVisitor";

  void serialize_bond_visitor_(std::ostream& ostr) const;

 private:
  double energy_;
};

inline std::shared_ptr<BondVisitor> MakeBondVisitor(
    const argtype &args = argtype()) {
  return std::make_shared<BondVisitor>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_BOND_VISITOR_H_
