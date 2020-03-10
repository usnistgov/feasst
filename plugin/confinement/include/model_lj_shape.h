
#ifndef FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_

#include "system/include/model_one_body.h"
#include "confinement/include/shape.h"
#include "math/include/constants.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  Note that the input shape of this model represents the shape of the cavity.

  /f$ U(r) = \epsilon \left( \frac{r}{\sigma} \right)^\alpha /f$
 */
class ModelLJShape : public ModelOneBody,
                       public ShapedEntity {
 public:
  ModelLJShape() {} // serialization only

  // Constructor
  ModelLJShape(
    /**
      alpha: set the exponent (default: 3).
     */
    std::shared_ptr<Shape> shape,
    const argtype& args = argtype())
    : ModelOneBody(), ShapedEntity(shape) {
    args_.init(args);
    alpha_ = args_.key("alpha").dflt("3").dble();
  }

  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const override {
    const int type = site.type();
    const double sigma = model_params.sigma().value(type);
    const double epsilon = model_params.epsilon().value(type);
    const double distance = shape()->nearest_distance(site.position());
    return epsilon * pow(distance/sigma, alpha_);
  }

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJShape>(istr); }
  explicit ModelLJShape(std::istream& istr);
  virtual ~ModelLJShape() {}

 private:
  const std::string class_name_ = "ModelLJShape";
  double alpha_;
  Arguments args_;
};

inline std::shared_ptr<ModelLJShape> MakeModelLJShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelLJShape>(shape);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
