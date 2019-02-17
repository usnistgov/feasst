
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
  virtual void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);

  // If model parameters are not given, then obtain them from the configuration.
  void compute(
      const ModelOneBody& model,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->unique_types().model_params();
    compute(model, model_params, config, group_index);
  }
  void compute(
      const ModelTwoBody& model,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->unique_types().model_params();
    compute(model, model_params, config, group_index);
  }
  void compute(
      const ModelTwoBody& model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->unique_types().model_params();
    compute(model, model_params, selection, config, group_index);
  }
  void compute(
      const ModelOneBody& model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->unique_types().model_params();
    compute(model, model_params, selection, config, group_index);
  }

  /// Return the energy.
  double energy() const { return energy_; }

  /// Set the energy.
  void set_energy(const double energy) { energy_ = energy; }

  virtual ~VisitModel() {}

  /// Test if energy of whole system is consistent with sum of energy
  /// of selection by particles.
  void check_energy(
      const Model& model,
      Configuration * config,
      const int group_index = 0);

  virtual void revert() {}

 protected:
  void inner_(
    const Site& site1,
    const Site& site2,
    const Domain& domain,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    Position * relative);

 private:
  double energy_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
