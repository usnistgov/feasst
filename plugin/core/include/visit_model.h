
#ifndef FEASST_CORE_VISIT_MODEL_H_
#define FEASST_CORE_VISIT_MODEL_H_

#include "core/include/model.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody;
class ModelTwoBody;
class ModelThreeBody;

class VisitModelInner {
 public:
  virtual void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    Position * relative);

  virtual void precompute(Configuration * config) {}

  void set_energy(const double energy) { energy_ = energy; }
  void add_energy(const double energy) { energy_ += energy; }
  double energy() const { return energy_; }

  virtual ~VisitModelInner() {}

 private:
  double energy_ = 0.;
};

/**
  See Model for a description of the compute methods. These are mirrored by
  simply switching the calling object and the first argument
  (.e.g, Model.compute(Visitor, ...) vs Visitor.compute(Model, ...)
 */
class VisitModel {
 public:
  VisitModel() {
    set_inner();
  }

  void set_inner(const std::shared_ptr<VisitModelInner> inner =
    std::make_shared<VisitModelInner>()) {
    inner_ = inner; }

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

  void zero_energy() {
    energy_ = 0.;
    inner_->set_energy(0.);
  }

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

  virtual void precompute(Configuration * config) {
    inner_->precompute(config); }

  const std::shared_ptr<VisitModelInner> inner() const { return inner_; }

 private:
  double energy_;
  std::shared_ptr<VisitModelInner> inner_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
