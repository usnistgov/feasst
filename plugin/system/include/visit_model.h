
#ifndef FEASST_SYSTEM_VISIT_MODEL_H_
#define FEASST_SYSTEM_VISIT_MODEL_H_

#include "system/include/model.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_inner.h"

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
  VisitModel() {
    set_inner(); }

  VisitModel(std::shared_ptr<VisitModelInner> inner) {
    set_inner(inner); }

  void set_inner(const std::shared_ptr<VisitModelInner> inner =
    std::make_shared<VisitModelInner>()) {
    inner_ = inner; }

  const VisitModelInner * const inner() const { return inner_.get(); }

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

  /// Test if energy of whole system is consistent with sum of energy
  /// of selection by particles.
  void check_energy(
      const Model& model,
      Configuration * config,
      const int group_index = 0);

  virtual void prep_for_revert() { inner_->prep_for_revert(); }
  virtual void revert() { inner_->revert(); }
  virtual void finalize() { inner_->finalize(); }

  virtual void precompute(Configuration * config) {
    inner_->precompute(config); }

  // serialization
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const {
    ostr << class_name_ << " ";
    serialize_visit_model_(ostr);
  }

  virtual std::shared_ptr<VisitModel> create(std::istream& istr) const {
    return std::make_shared<VisitModel>(istr);
  }

  std::map<std::string, std::shared_ptr<VisitModel> >& deserialize_map();

  std::shared_ptr<VisitModel> deserialize(std::istream& istr) {
    return template_deserialize(deserialize_map(), istr); }

  VisitModel(std::istream& istr);
  virtual ~VisitModel() {}

 protected:
  void serialize_visit_model_(std::ostream& ostr) const;
//  void deserialize_visit_model_(std::istream& istr, std::shared_ptr<VisitModel> visitor) const;
  VisitModelInner * get_inner_() const { return inner_.get(); }

  // HWH hacky addition: also, prep inner for reverting,
  // because this is called at beginning of every pair-wise selection compute
  // optimization to avoid repeated construction of Position.
  Position relative_;
  void init_relative_(const Domain& domain, Position * relative) {
    prep_for_revert();
    if (relative->dimension() != domain.dimension()) {
      relative->set_vector(domain.side_length().coord());
    }
  }

 private:
  const std::string class_name_ = "VisitModel";
  double energy_ = 0.;
  std::shared_ptr<VisitModelInner> inner_;
};

inline std::shared_ptr<VisitModel> MakeVisitModel() {
  return std::make_shared<VisitModel>();
}

inline std::shared_ptr<VisitModel> MakeVisitModel(
    std::shared_ptr<VisitModelInner> inner) {
  return std::make_shared<VisitModel>(inner);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_H_
