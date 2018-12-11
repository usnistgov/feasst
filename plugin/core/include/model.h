
#ifndef FEASST_CORE_MODEL_H_
#define FEASST_CORE_MODEL_H_

#include "core/include/particle.h"
#include "core/include/configuration.h"

namespace feasst {

class VisitModel;

class Model {
 public:
  /// HWH depreciate
  virtual double compute(VisitModel& visitor, const Configuration& config, const int iPart) = 0;

  virtual double compute_selection(VisitModel& visitor,
                                   const Configuration& config) = 0;
  virtual ~Model() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_H_
