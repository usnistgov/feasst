
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
  VisitModelInner() {}

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

  // serialize
  virtual void serialize(std::ostream& ostr) const {
    serialize_visit_model_inner_(ostr); }

  VisitModelInner(std::istream& istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize(&energy_, istr);
  }

  virtual std::shared_ptr<VisitModelInner> create(std::istream& istr) const {
    return std::make_shared<VisitModelInner>(istr);
  }

  std::map<std::string, std::shared_ptr<VisitModelInner> >& deserialize_map();

  std::shared_ptr<VisitModelInner> deserialize(std::istream& istr) {
    return template_deserialize(deserialize_map(), istr); }

  virtual ~VisitModelInner() {}

 protected:
  void serialize_visit_model_inner_(std::ostream& ostr) const {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize(energy_, ostr);
  }


 private:
  const std::string class_name_ = "VisitModelInner";
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

  // serialization
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

 private:
  const std::string class_name_ = "VisitModel";
  double energy_ = 0.;
  std::shared_ptr<VisitModelInner> inner_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
