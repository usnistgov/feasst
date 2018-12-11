
#ifndef FEASST_CORE_MODEL_ONE_BODY_H_
#define FEASST_CORE_MODEL_ONE_BODY_H_

#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/particle.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody : public Model {
 public:
  // HWH depreciate
  double compute(VisitModel& visitor, const Configuration& config, const int iPart) override {
    visitor.loop_by_particle(config, *this, iPart);
    return visitor.energy();
  }
  double compute_selection(VisitModel& visitor,
                           const Configuration& config) override {
    visitor.energy_of_selection(config, *this);
    return visitor.energy();
  }
  virtual double evaluate(const Site& site,
                          const Configuration& config,
                          const ModelParams& model_params) const = 0;
  virtual ~ModelOneBody() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_ONE_BODY_H_
