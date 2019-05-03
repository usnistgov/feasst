
#ifndef FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_

#include "core/include/model_one_body.h"
#include "confinement/include/shape.h"
#include "core/include/constants.h"

namespace feasst {

/**
  Note that the input shape of this model represents the shape of the cavity.
  Thus, the hard interaction is outside of the shape, not inside.
 */
class ModelHardShape : public ModelOneBody,
                       public ShapedEntity {
 public:
  ModelHardShape() {} // serialization only
  ModelHardShape(std::shared_ptr<Shape> shape)
    : ModelOneBody(), ShapedEntity(shape) {}

  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const override {
    const int type = site.type();
    const double sigma = model_params.sigma().value(type);
    if (shape()->is_inside(site.position(), sigma)) {
      return 0.;
    }
    return NEAR_INFINITY;
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    ShapedEntity::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelHardShape>(istr);
  }

  ModelHardShape(std::istream& istr) : ModelOneBody(), ShapedEntity(istr) {
    feasst_deserialize_version(istr);
  }

  virtual ~ModelHardShape() {}

 private:
  const std::string class_name_ = "ModelHardShape";
};

inline std::shared_ptr<ModelHardShape> MakeModelHardShape(
    std::shared_ptr<Shape> shape) {
  return std::make_shared<ModelHardShape>(shape);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
