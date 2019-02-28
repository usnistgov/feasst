
#ifndef FEASST_CORE_POTENTIAL_H_
#define FEASST_CORE_POTENTIAL_H_

#include <vector>
#include <memory>
#include <string>
#include "core/include/debug.h"
#include "core/include/configuration.h"
#include "core/include/visit_model.h"
#include "core/include/model_empty.h"
#include "core/include/constants.h"

namespace feasst {

/**
  A potential represents both the model and the method used to compute the model
  for a given configuration.
  A potential also allows customization of the model parameters.
 */
class Potential {
 public:
  Potential() {
    model_ = std::make_shared<ModelEmpty>();
    visit_model_ = std::make_shared<VisitModel>();
    set_group_index(); }

  /// Set the index of the group of configuration which contributes to this
  /// potential. By default, the index is zero, which represents the entire
  /// configuration.
  /// Note that for VisitModelCell, group_index refers to the cell_index.
  void set_group_index(const int group_index = 0) {group_index_ = group_index; }

  /// Return the index of the group
  int group_index() const { return group_index_; }

  /// Set the method used to compute.
  void set_visit_model(std::shared_ptr<VisitModel> visit_model) {
    visit_model_ = visit_model; }

  /// Return the method used to compute.
  const std::shared_ptr<VisitModel> visit_model() const {
    return visit_model_; }

  /// Set the model.
  void set_model(std::shared_ptr<Model> model) { model_ = model; }

  /// Return the model.
  const std::shared_ptr<Model> model() const { return model_; }

  /// Set the model parameters. If not set, use the one from configuration.
  void set_model_params(const ModelParams& model_params) {
    model_params_override_ = true;
    model_params_ = model_params;
  }

  /// Set the model parameters to the one in the configuration.
  void set_model_params(const Configuration& config) {
    set_model_params(config.model_params());
  }

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const char* name,
                       const int site_type,
                       const double value) {
    ASSERT(model_params_override_, "you must first initialize model params");
    model_params_.set(name, site_type, value);
  }

  /// Return the model parameters.
  const ModelParams& model_params() const {
    ASSERT(model_params_override_, "you must first initialize model params");
    return model_params_;
  }

  /// Compute the energy of the entire configuration.
  double energy(Configuration * config) {
    ASSERT(visit_model_, "visitor must be set.");
    if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, group_index_, config, visit_model_.get());
    } else {
      stored_energy_ = model_->compute(group_index_, config, visit_model_.get());
    }
    return stored_energy_;
  }

  /// Compute the energy of a selection of the configuration.
  double energy(const Select& select, Configuration * config) {
    ASSERT(visit_model_, "visitor must be set.");
    if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, select, group_index_, config, visit_model_.get());
    } else {
      stored_energy_ = model_->compute(select, group_index_, config, visit_model_.get());
    }
    return stored_energy_;
  }

  /// Return the last computed value of the energy.
  double stored_energy() const { return stored_energy_; }

  // Revert any changes to the configuration due to the last energy computation.
  void revert() { visit_model_->revert(); }

 private:
  int group_index_;
  std::shared_ptr<VisitModel> visit_model_;
  std::shared_ptr<Model> model_;
  double stored_energy_;

  bool model_params_override_ = false;
  ModelParams model_params_;
};

}  // namespace feasst

#endif  // FEASST_CORE_POTENTIAL_H_
