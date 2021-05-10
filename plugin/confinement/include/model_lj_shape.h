
#ifndef FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_

#include "system/include/model_one_body.h"
#include "shape/include/shape.h"
#include "utils/include/arguments.h"

namespace feasst {

class ModelLJShapeEnergyAtCutoff;

/**
  Note that the input shape of this model represents the shape of the cavity.

  \f$
  U_{LJ}(r) = \epsilon \left( \frac{\sigma}{r} \right)^\alpha
  U_{LJ}^{CS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) & : r < r_c
      0 & : r \ge r_c
    \end{array}
  \right. \f$

  where \f$r\f$ is the nearest distance of the site to the confinement.
  Note that the potential is cut and shifted to zero at the cut off distance,
  \f$r_c\f$.
 */
class ModelLJShape : public ModelOneBody,
                     public ShapedEntity {
 public:
  // for serialization only
  ModelLJShape() { class_name_ = "ModelLJShape"; }

  /**
    args:
    - alpha: set the exponent (default: 3).
    - disable_shift: disable shifting of the potential to zero (default: false).
   */
  ModelLJShape(
    std::shared_ptr<Shape> shape,
    argtype args = argtype());

  ModelLJShape(std::shared_ptr<Shape> shape,
    argtype * args);

  /// Precompute the shift factor for optimization.
  void precompute(const ModelParams& existing) override;

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  double energy(const double epsilon, const double sigma,
    const double distance) const;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJShape>(istr); }
  explicit ModelLJShape(std::istream& istr);
  virtual ~ModelLJShape() {}

 private:
  double alpha_;
  bool disable_shift_;
  std::shared_ptr<ModelLJShapeEnergyAtCutoff> shift_;
};

inline std::shared_ptr<ModelLJShape> MakeModelLJShape(
    std::shared_ptr<Shape> shape, argtype args = argtype()) {
  return std::make_shared<ModelLJShape>(shape, args);
}

class ModelLJShapeEnergyAtCutoff : public ModelParam {
 public:
  ModelLJShapeEnergyAtCutoff() : ModelParam() {}

  void set_model(ModelLJShape * model) { model_ = model; }

  double compute(const int type1, const ModelParams& model_params) override {
    const double epsilon = model_params.epsilon().value(type1);
    const double sigma = model_params.sigma().value(type1);
    const double cutoff = model_params.cutoff().value(type1);
    return model_->energy(epsilon, sigma, cutoff);
  }

  ModelLJShapeEnergyAtCutoff(std::istream& istr) : ModelParam(istr) {}

private:
  /// temporary pointer for use with compute without requiring new arguments.
  ModelLJShape * model_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
