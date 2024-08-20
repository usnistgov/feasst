
#ifndef FEASST_SYSTEM_VISIT_MODEL_BOND_H_
#define FEASST_SYSTEM_VISIT_MODEL_BOND_H_

#include <map>
#include <string>
#include <memory>
#include "system/include/visit_model.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Interactions between bonded sites using "inter"-like models are computed.
 */
class VisitModelBond : public VisitModel {
 public:
  explicit VisitModelBond(argtype args = argtype());
  explicit VisitModelBond(argtype * args);
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<VisitModelBond>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<VisitModelBond>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelBond(std::istream& istr);
  virtual ~VisitModelBond() {}
};

inline std::shared_ptr<VisitModelBond> MakeVisitModelBond(
    const argtype& args = argtype()) {
  return std::make_shared<VisitModelBond>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_BOND_H_
