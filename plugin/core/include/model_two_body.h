
#ifndef FEASST_CORE_MODEL_TWO_BODY_H_
#define FEASST_CORE_MODEL_TWO_BODY_H_

#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/particle.h"

namespace feasst {

class ModelTwoBody : public Model {
 public:
  double compute(VisitModel& visitor,
      const Configuration& config,
      const int group_index) const override {
    visitor.compute(config, *this, group_index);
    return visitor.energy();
  }
  double compute(VisitModel& visitor,
      const Configuration& config) const override {
    visitor.compute(config, *this, 0);
    return visitor.energy();
  }
  double compute(VisitModel& visitor,
      const Configuration& config,
      const Select& selection,
      const int group_index) const override {
    visitor.compute(config, *this, selection, group_index);
    return visitor.energy();
  }
  double compute(VisitModel& visitor,
      const Configuration& config,
      const Select& selection) const override {
    visitor.compute(config, *this, selection, 0);
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
