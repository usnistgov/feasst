
#ifndef FEASST_MODEL_EXPANDED_MODEL_EXPANDED_H_
#define FEASST_MODEL_EXPANDED_MODEL_EXPANDED_H_

#include <memory>
#include <string>
#include <vector>
#include "system/include/model_two_body_factory.h"

namespace feasst {

/**
  Contains a collection of two body models.
  This allows the visitor to consider multiple two body models for each
  pairwise interaction.
 */
class ModelExpanded : public ModelTwoBodyFactory {
 public:
  //@{
  /** @name Arguments
    - model_index: index of current model (default: 0).
    - ModelTwoBodyFactory arguments.
   */
  explicit ModelExpanded(argtype args = argtype());
  explicit ModelExpanded(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the index of the model.
  int model_index() const override { return model_index_; }

  /// Set the index of the model.
  void set_model_index(const int index) override;

  /// Return the energy of only the model with the particular index.
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  // serialize
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelExpanded>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelExpanded>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelExpanded(std::istream& istr);

  //@}
 private:
  int model_index_;
};

inline std::shared_ptr<ModelExpanded> MakeModelExpanded(argtype args = argtype()) {
  return std::make_shared<ModelExpanded>(args);
}

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_MODEL_EXPANDED_H_
