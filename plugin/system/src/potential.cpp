#include <memory>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/potential.h"
#include "system/include/model_empty.h"

namespace feasst {

Potential::Potential(const argtype& args) {
  model_ = std::make_shared<ModelEmpty>();
  visit_model_ = std::make_shared<VisitModel>();

  Arguments args_(args);
  group_index_ = args_.key("group_index").dflt("0").integer();
  if (args_.key("cell_index").used()) {
    ASSERT(group_index_ == 0, "cell_index overrides group_index");
    group_index_ = args_.integer();
  }
  prevent_cache_ = args_.key("prevent_cache").dflt("False").boolean();
}

Potential::Potential(std::shared_ptr<Model> model,
                     const argtype& args) : Potential(args) {
  model_ = model;
}

Potential::Potential(std::shared_ptr<VisitModel> visit_model,
                     const argtype& args) : Potential(args) {
  visit_model_ = visit_model;
}

Potential::Potential(
    std::shared_ptr<Model> model,
    std::shared_ptr<VisitModel> visit_model,
    const argtype& args) : Potential(model, args) {
  visit_model_ = visit_model;
}

void Potential::set(const ModelParams& model_params) {
  model_params_override_ = true;
  model_params_ = ModelParams(model_params);
}

void Potential::set_model_param(const char* name,
    const int site_type,
    const double value) {
  ASSERT(model_params_override_, "you must first initialize model params "
    << "before setting them.");
  model_params_.set(name, site_type, value);
}

const ModelParams& Potential::model_params() const {
  ASSERT(model_params_override_, "When model parameters are not overriden, "
    << "you must also provide the configuration as an argument.");
  return model_params_;
}

const ModelParams& Potential::model_params(const Configuration& config) const {
  if (model_params_override_) {
    return model_params_;
  }
  return config.model_params();
}

double Potential::energy(Configuration * config) {
  ASSERT(visit_model_, "visitor must be set.");
  if (prevent_cache_ || !cache_.is_unloading(&stored_energy_)) {
    if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, group_index_, config,
                                       visit_model_.get());
    } else {
      stored_energy_ = model_->compute(group_index_, config, visit_model_.get());
    }
    cache_.load(stored_energy_);
  }
  return stored_energy_;
}

double Potential::energy(const Select& select, Configuration * config) {
  ASSERT(visit_model_, "visitor must be set.");
  if (prevent_cache_ || !cache_.is_unloading(&stored_energy_)) {
    if (model_params_override_) {
      stored_energy_ = model_->compute(model_params_, select, group_index_,
                                       config, visit_model_.get());
    } else {
      stored_energy_ = model_->compute(select, group_index_, config,
                                       visit_model_.get());
    }
    cache_.load(stored_energy_);
  }
  return stored_energy_;
}

int Potential::cell_index() const {
  ASSERT(visit_model_->class_name() == "VisitModelCell", "error");
  return group_index();
}

void Potential::precompute(Configuration * config) {
  visit_model_->precompute(config);
  const ModelParams& params = model_params(*config);
  model_->precompute(params);
  const double max_cutoff = maximum(params.cutoff().values());
  const double half_min_side = 0.5*config->domain().min_side_length();
  if (max_cutoff - NEAR_ZERO > half_min_side) {
    WARN("The maximum cutoff:" << max_cutoff << " is greater than half the " <<
         "minimum side length: " << half_min_side);
  }
}

void Potential::check() const {
  visit_model_->check();
}

void Potential::serialize(std::ostream& ostr) const {
  feasst_serialize_version(432, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstdr(visit_model_, ostr);
  feasst_serialize_fstdr(model_, ostr);
  feasst_serialize(stored_energy_, ostr);
  feasst_serialize(model_params_override_, ostr);
  if (model_params_override_) {
    feasst_serialize_fstobj(model_params_, ostr);
  }
  feasst_serialize_fstobj(cache_, ostr);
  feasst_serialize(prevent_cache_, ostr);
}

Potential::Potential(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(432 == version, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  // feasst_deserialize_fstdr(visit_model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      visit_model_ = visit_model_->deserialize(istr);
    }
  }
  // feasst_deserialize_fstdr(model_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      model_ = model_->deserialize(istr);
    }
  }
  feasst_deserialize(&stored_energy_, istr);
  feasst_deserialize(&model_params_override_, istr);
  if (model_params_override_) {
    feasst_deserialize_fstobj(&model_params_, istr);
  }
  feasst_deserialize_fstobj(&cache_, istr);
  feasst_deserialize(&prevent_cache_, istr);
}

void Potential::set_model_params(const Configuration& config) {
  set(config.model_params());
}

}  // namespace feasst
