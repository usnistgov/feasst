
#ifndef FEASST_CORE_VISIT_MODEL_H_
#define FEASST_CORE_VISIT_MODEL_H_

#include "core/include/model.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody;
class ModelTwoBody;
class ModelThreeBody;

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
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
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
      const ModelThreeBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) {
    ERROR("not implemented");
  }

  // If model parameters are not given, then obtain them from the configuration.
  void compute(
      const ModelOneBody& model,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->model_params();
    compute(model, model_params, config, group_index);
  }
  void compute(
      const ModelOneBody& model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->model_params();
    compute(model, model_params, selection, config, group_index);
  }
  void compute(
      const ModelTwoBody& model,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->model_params();
    compute(model, model_params, config, group_index);
  }
  void compute(
      const ModelTwoBody& model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->model_params();
    compute(model, model_params, selection, config, group_index);
  }
  void compute(
      const ModelThreeBody& model,
      Configuration * config,
      const int group_index = 0) {
    const ModelParams& model_params = config->model_params();
    compute(model, model_params, config, group_index);
  }

  /// Return the energy.
  double energy() const { return energy_; }

  /// Set the energy.
  void set_energy(const double energy) { energy_ = energy; }

  /// Increment the energy.
  void increment_energy(const double energy) { energy_ += energy; }

  virtual ~VisitModel() {}

  /// Test if energy of whole system is consistent with sum of energy
  /// of selection by particles.
  void check_energy(
      const Model& model,
      Configuration * config,
      const int group_index = 0);

  virtual void revert() {}

 protected:
  virtual void inner_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    Position * relative);

 private:
  double energy_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
