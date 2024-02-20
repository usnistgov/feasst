
#ifndef FEASST_MODEL_EXPANDED_PERTURB_MODEL_H_
#define FEASST_MODEL_EXPANDED_PERTURB_MODEL_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  Change the ModelExpanded::model_index by plus or minus one.
 */
class PerturbModel : public Perturb {
 public:
  //@{
  /** @name Arguments
    args:
    - potential_index: index of potential with model_index (default: 0).
   */

  /**
   */
  explicit PerturbModel(argtype args = argtype());
  explicit PerturbModel(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held = false,
    Acceptance * acceptance = NULL) override;

  /// Change Model::model_index.
  void change_model_index(const int delta_model_index, System * system);

  void revert(System * system) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbModel(std::istream& istr);
  virtual ~PerturbModel() {}
  //@}

 private:
  int potential_index_;

  // temporary and not serialized
  int previous_model_;
};

inline std::shared_ptr<PerturbModel> MakePerturbModel(argtype args = argtype()) {
  return std::make_shared<PerturbModel>(args);
}

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_PERTURB_MODEL_H_
