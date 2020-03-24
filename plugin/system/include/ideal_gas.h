
#ifndef FEASST_SYSTEM_IDEAL_GAS_H_
#define FEASST_SYSTEM_IDEAL_GAS_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The ideal gas does not interact.
  Thus, this model returns 0 in call cases.
 */
class IdealGas : public ModelTwoBody {
 public:
  IdealGas() {}

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    return 0.;
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<IdealGas>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit IdealGas(std::istream& istr);
  virtual ~IdealGas() {}

 protected:
  void serialize_ideal_gas_(std::ostream& ostr) const;

 private:
  const std::string class_name_ = "IdealGas";
};

inline std::shared_ptr<IdealGas> MakeIdealGas() {
  return std::make_shared<IdealGas>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_IDEAL_GAS_H_
