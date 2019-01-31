
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
  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const {
    const int type = site.type();
    const double sigma = model_params.sigma().value(type);
    if (shape()->is_inside(site.position(), sigma)) {
      return 0.;
    }
    return NEAR_INFINITY;
  }

  virtual ~ModelHardShape() {}
 private:
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_HARD_SHAPE_H_
