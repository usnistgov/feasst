
#ifndef FEASST_CORE_VISIT_MODEL_H_
#define FEASST_CORE_VISIT_MODEL_H_

#include "core/include/model.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody;
class ModelTwoBody;

/**
  See Model for a description of the compute methods. These are mirrored by
  simply switching the calling object and the first argument
  (.e.g, Model.compute(Visitor, ...) vs Visitor.compute(Model, ...)
 */
class VisitModel {
 public:
  virtual void compute(const Configuration& config,
      const ModelOneBody& model,
      const int group_index = 0);
  virtual void compute(const Configuration& config,
      const ModelTwoBody& model,
      const int group_index = 0);
  virtual void compute(const Configuration& config,
      const ModelTwoBody& model,
      const Select& selection,
      const int group_index = 0);
  void compute(const Configuration& config,
      const ModelOneBody& model,
      const Select& selection,
      const int group_index = 0);

  /// Return the energy.
  double energy() const { return energy_; }

  /// Set the energy.
  void set_energy(const double energy) { energy_ = energy; }

  virtual ~VisitModel() {}

 private:
  double energy_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
