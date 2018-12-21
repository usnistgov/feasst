
#ifndef FEASST_CORE_SYSTEM_H_
#define FEASST_CORE_SYSTEM_H_

#include <vector>
#include "core/include/debug.h"
#include "core/include/configuration.h"
#include "core/include/visit_model.h"
#include "core/include/model_lj.h"

namespace feasst {
/**
  Systems may have multiple configurations but their typing and grouping should be the same.
  HWH refactor how the configurations are set up (e.g., no add_configuration).
  This way we can enforce typing.
  Allow duplication of configuration.
  Or maybe this should be done in the configuration class itself?
 */
class System {
 public:
  /// Set the configuration.
  void add_configuration(const Configuration& configuration) { configurations_.push_back(configuration); }

  /// Return the configuration
  const Configuration& configuration(const int index = 0) const { return configurations_[index]; }
  Configuration* get_configuration(const int index = 0) { return &configurations_[index]; }

  int dimension() const { return configurations_.front().dimension(); }

  double energy() {
    return energy(model_);
  }

  double energy(const Select& selection) {
    return model_.compute(visit_model_, configurations_.front(), selection);
  }

  double energy(const ModelTwoBody& model) {
    ASSERT(dimension() > 0, "size error");
    visit_model_.compute(configurations_.front(), model);
    return visit_model_.energy();
  }

  VisitModel* get_visitor() { return &visit_model_; }
  const VisitModel& visitor() const { return visit_model_; }

 private:
  std::vector<Configuration> configurations_;
  VisitModel visit_model_;
  ModelLJ model_;
  /// ModelLJSingleComp model_;
//  vector<ModelTwoBody> two_body_models_;
//  vector<ModelOneBody> one_body_models_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_H_
