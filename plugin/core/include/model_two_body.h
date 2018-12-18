
#ifndef FEASST_CORE_MODEL_TWO_BODY_H_
#define FEASST_CORE_MODEL_TWO_BODY_H_

#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/particle.h"

namespace feasst {

class ModelTwoBody : public Model {
 public:
  /// HWH depreciate
  double compute(
      VisitModel& visitor,
      const Configuration& config,
      const int iPart) override {
    visitor.loop_by_particle(config, *this, iPart);
    return visitor.energy();
  }
  double compute(
      VisitModel& visitor,
      const Configuration& config,
      const Select& selection) override {
    visitor.compute(config, *this, selection);
    return visitor.energy();
  }
  virtual double evaluate(
      const Position &relative,
      const Site& site1,
      const Site& site2,
      const ModelParams& model_params) const = 0;
  virtual ~ModelTwoBody() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_TWO_BODY_H_
