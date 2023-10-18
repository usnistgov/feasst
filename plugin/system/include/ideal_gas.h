
#ifndef FEASST_SYSTEM_IDEAL_GAS_H_
#define FEASST_SYSTEM_IDEAL_GAS_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The ideal gas does not interact.
  Thus, this model returns 0 energy in call cases.
 */
class IdealGas : public ModelTwoBody {
 public:
  IdealGas() { class_name_ = "IdealGas"; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override {
    return 0.;
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<IdealGas>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<IdealGas>(); }
  void serialize(std::ostream& ostr) const override;
  explicit IdealGas(std::istream& istr);
  virtual ~IdealGas() {}

 protected:
  void serialize_ideal_gas_(std::ostream& ostr) const;
};

inline std::shared_ptr<IdealGas> MakeIdealGas() {
  return std::make_shared<IdealGas>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_IDEAL_GAS_H_
