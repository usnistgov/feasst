
#ifndef FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
#define FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_

#include "utils/include/arguments.h"
#include "configuration/include/model_params.h"
#include "system/include/model_one_body.h"
#include "shape/include/shape.h"

namespace feasst {

class ModelLJShapeEnergyAtCutoff;

/**
  Note that the input shape of this model represents the shape of the cavity.

  \f$U_{LJ}(r) = \epsilon \left( \frac{\sigma}{r + \Delta} \right)^\alpha\f$

  \f$U_{LJ}^{CS}(r) = \left\{
    \begin{array}{lr}
      U_{LJ}(r) - U_{LJ}(r_c) & : r < r_c \\
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

  //@{
  /** @name Arguments
    - Same arguments as ShapeFile.
    - alpha: set the exponent (default: 3.).
    - delta: set the delta parameter (default: 0.).
    - disable_shift: disable shifting of the potential to zero (default: false).
    - wall_sigma: If != 0 (default: 0), use Lorentz-Berthelot mixing rules
      between this wall sigma and the fluid sigma.
      Otherwise, the sigma for each site type may be set with the argument
      Potential::sigma[i].
    - wall_epsilon: If != 0 (default: 0), use a slighly modified version of
      Lorentz-Berthelot mixing rules between this wall epsilon and the fluid.
      To account for negative epsilon (attractions), the mixing rule is
      sign(wall_epsilon)sqrt(|wall_epsilon|*fluid_epsilon).
      Otherwise, the epsilon for each site type may be set with the argument
      Potential::epsilon[i].
   */
  explicit ModelLJShape(argtype args);
  explicit ModelLJShape(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  explicit ModelLJShape(std::shared_ptr<Shape> shape, argtype args = argtype());
  ModelLJShape(std::shared_ptr<Shape> shape, argtype * args);

  /// Precompute the shift factor for optimization.
  void precompute(const ModelParams& existing) override;

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  /// Return the sigma (in case optionally mixed)
  double sigma(const int site_type, const ModelParams& params);

  /// Return the epsilon (in case optionally mixed)
  double epsilon(const int site_type, const ModelParams& params);

  double energy(const double epsilon, const double sigma,
    const double distance) const;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelLJShape>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelLJShape>(args); }
  explicit ModelLJShape(std::istream& istr);
  virtual ~ModelLJShape() {}

  //@}
 private:
  double alpha_;
  double delta_;
  bool disable_shift_;
  std::shared_ptr<ModelLJShapeEnergyAtCutoff> shift_;
  double wall_sigma_;
  Sigma mixed_sigma_;
  double wall_epsilon_;
  Epsilon mixed_epsilon_;

  void parse_args_(argtype * args);
};

inline std::shared_ptr<ModelLJShape> MakeModelLJShape(
    std::shared_ptr<Shape> shape, argtype args = argtype()) {
  return std::make_shared<ModelLJShape>(shape, args);
}

class ModelLJShapeEnergyAtCutoff : public ModelParam {
 public:
  ModelLJShapeEnergyAtCutoff() : ModelParam() {}
  void set_model(ModelLJShape * model) { model_ = model; }
  double compute(const int type1, const ModelParams& model_params) override;
  ModelLJShapeEnergyAtCutoff(std::istream& istr) : ModelParam(istr) {}

private:
  /// temporary pointer for use with compute without requiring new arguments.
  ModelLJShape * model_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_LJ_SHAPE_H_
