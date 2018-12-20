
#ifndef FEASST_CORE_MODEL_ONE_BODY_H_
#define FEASST_CORE_MODEL_ONE_BODY_H_

#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/particle.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody : public Model {
 public:
  double compute(VisitModel& visitor,
      const Configuration& config,
      const int group_index) override {
    visitor.compute(config, *this, group_index);
    return visitor.energy();
  }
  double compute(VisitModel& visitor, const Configuration& config) override {
    visitor.compute(config, *this, 0);
    return visitor.energy();
  }
  double compute(VisitModel& visitor,
      const Configuration& config,
      const Select& selection,
      const int group_index) override {
    visitor.compute(config, *this, selection, group_index);
    return visitor.energy();
  }
  double compute(VisitModel& visitor,
      const Configuration& config,
      const Select& selection) override {
    visitor.compute(config, *this, selection, 0);
    return visitor.energy();
  }
  virtual double evaluate(
      const Site& site,
      const Configuration& config,
      const ModelParams& model_params) const = 0;
  virtual ~ModelOneBody() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_ONE_BODY_H_
