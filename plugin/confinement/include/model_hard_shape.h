
#ifndef FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_

#include "system/include/model_one_body.h"
#include "shape/include/shape.h"

namespace feasst {

/**
  Create a cavity of a given shape, where the potential energy outside of the
  shape is infinite.
  Optionally set the cavity to false (default: true) to model a hard shape.
 */
class ModelHardShape : public ModelOneBody,
                       public ShapedEntity {
 public:
  // for serialization only
  ModelHardShape() { class_name_ = "ModelHardShape"; }

  /**
    args:
    - Same arguments as ShapeFile.
    - cavity: if true, the potential outside the shape is infinite.
      if false, the potential inside the shape is infinite (default: true).
   */
  explicit ModelHardShape(argtype args);
  explicit ModelHardShape(argtype * args);

  ModelHardShape(std::shared_ptr<Shape> shape)
    : ModelOneBody(), ShapedEntity(shape) { class_name_ = "ModelHardShape"; }

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelHardShape>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelHardShape>(args); }
  explicit ModelHardShape(std::istream& istr);
  virtual ~ModelHardShape() {}

 private:
  bool cavity_ = true;
};

inline std::shared_ptr<ModelHardShape> MakeModelHardShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelHardShape>(shape);
}

inline std::shared_ptr<ModelHardShape> MakeModelHardShape(argtype args) {
  return std::make_shared<ModelHardShape>(args); }

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
