
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_

#include <memory>
#include <string>
#include <vector>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  Contains a collection of two body models.
  This allows the visitor to consider multiple two body models for each
  pairwise interaction.
 */
class ModelTwoBodyFactory : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - model[i]: add the i-th model.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - model_file: if not empty, file name that lists the ModelTwoBody
      (default: empty).
      The first line is simply ModelTwoBodyFactory.
      The second line is empty.
      The remaining lines are for creating each ModelTwoBody.
   */
  explicit ModelTwoBodyFactory(argtype args = argtype());
  explicit ModelTwoBodyFactory(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Add a vector of two body models.
  void add(std::vector<std::shared_ptr<ModelTwoBody> > models);

  /// Add a two body model.
  void add(std::shared_ptr<ModelTwoBody> model) {
    models_.push_back(model); }

  /// Return the number of two body models.
  int num() const { return static_cast<int>(models_.size()); }

  /// Return a two body model by index.
  const Model * model(const int index) const { return models_[index].get(); }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const ModelParams& existing) override;

  // serialize
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTwoBodyFactory>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelTwoBodyFactory>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelTwoBodyFactory(std::istream& istr);

  //@}
 protected:
  void serialize_model_two_body_factory_(std::ostream& ostr) const;
  std::vector<std::shared_ptr<Model> > models_;
};

inline std::shared_ptr<ModelTwoBodyFactory> MakeModelTwoBodyFactory(argtype args = argtype()) {
  return std::make_shared<ModelTwoBodyFactory>(args);
}

std::shared_ptr<ModelTwoBodyFactory> MakeModelTwoBodyFactory(
    std::shared_ptr<ModelTwoBody> model1,
    std::shared_ptr<ModelTwoBody> model2);

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_
