#ifndef FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_
#define FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Attempt to change the Model::model_index by plus or minus one.
 */
class TrialModel : public Trial {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - PerturbModel arguments.
    - Tunable arguments.
   */
  explicit TrialModel(argtype args = argtype());
  explicit TrialModel(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialModel>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialModel>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialModel(std::istream& istr);
  virtual ~TrialModel() {}
  //@}
};

inline std::shared_ptr<TrialModel> MakeTrialModel(argtype args = argtype()) {
  return std::make_shared<TrialModel>(args); }

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_TRIAL_MODEL_H_
