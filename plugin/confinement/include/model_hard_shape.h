
#ifndef FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_

#include "system/include/model_one_body.h"
#include "shape/include/shape.h"

namespace feasst {

/**
  Create a cavity of a given shape, where the potential energy outside of the
  shape is infinite.
 */
class ModelHardShape : public ModelOneBody,
                       public ShapedEntity {
 public:
  // for serialization only
  ModelHardShape() { class_name_ = "ModelHardShape"; }
  ModelHardShape(std::shared_ptr<Shape> shape)
    : ModelOneBody(), ShapedEntity(shape) { class_name_ = "ModelHardShape"; }

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelHardShape>(istr); }
  explicit ModelHardShape(std::istream& istr);
  virtual ~ModelHardShape() {}
};

inline std::shared_ptr<ModelHardShape> MakeModelHardShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelHardShape>(shape);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
