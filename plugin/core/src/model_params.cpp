
#include <algorithm>
#include "core/include/model_params.h"
#include "core/include/debug.h"
#include "core/include/utils_math.h"

namespace feasst {

void ModelParam::add(const double value) {
  values_.push_back(value);
  max_value_ = *std::max_element(values_.begin(), values_.end());
}

void ModelParam::add(const Site site, const double default_value) {
  double value;
  if (site.properties().value(name(), &value)) {
    add(value);
  } else {
    add(default_value);
  }
}

void ModelParam::add(const Particle particle) {
  for (const Site& site : particle.sites()) {
    add(site);
  }
}

void ModelParam::mix() {
  ASSERT(!is_mixed_override_, "mixed values were set specifically. " <<
    "And they would be overridden by this call to mix. Set mixed values last");
  resize(size(), size(), &mixed_values_);
  for (int index1 = 0; index1 < size(); ++index1) {
    const double value1 = values_[index1];
    for (int index2 = 0; index2 < size(); ++index2) {
      const double value2 = values_[index2];
      mixed_values_[index1][index2] = mix_(value1, value2);
    }
  }
  max_mixed_value_ = maximum(mixed_values_);
}

double ModelParam::mixed_max() const {
  ASSERT(mixed_values_.size() > 0, "no max");
  return max_mixed_value_;
}

double ModelParam::value(const int type) const {
  ASSERT(type < size(),
    "size error: type(" << type << ") size(" << size() << ")");
  return values_[type];
}

double ModelParam::compute(const int type1, const int type2,
    const ModelParams& model_params) {
  ERROR("not implemented");
  return 0.;
}

void ModelParam::set_param(const ModelParams& existing) {
  // initialize the number of types with zero value
  for (int type = 0; type < existing.size(); ++type) {
    add(0.);
  }
  mix();

  // set the values of the mixed types
  for (int type1 = 0; type1 < existing.size(); ++type1) {
    for (int type2 = 0; type2 < existing.size(); ++type2) {
      set_mixed(type1, type2, compute(type1, type2, existing));
    }
  }
}

void ModelParams::add_() {
  add(epsilon_);
  add(sigma_);
  add(cutoff_);
  add(charge_);
}

ModelParams::ModelParams() {
  epsilon_ = std::make_shared<Epsilon>();
  sigma_ = std::make_shared<Sigma>();
  cutoff_ = std::make_shared<CutOff>();
  charge_ = std::make_shared<Charge>();
  add_();
}

ModelParams::ModelParams(const ModelParams& params) {
  epsilon_ = std::make_shared<Epsilon>(*params.epsilon_);
  sigma_ = std::make_shared<Sigma>(*params.sigma_);
  cutoff_ = std::make_shared<CutOff>(*params.cutoff_);
  charge_ = std::make_shared<Charge>(*params.charge_);
  add_();
  for (int extra = 4; extra < static_cast<int>(params.params_.size());
       ++extra) {
    add(params.params_[extra]);
  }
  set_properties(params.properties());
}

void ModelParams::add(const Particle particle) {
  for (std::shared_ptr<ModelParam> param : params_) {
    param->add(particle);
  }
  mix();
}

void ModelParams::mix() {
  for (std::shared_ptr<ModelParam> param : params_) {
    param->mix();
  }
}

int ModelParams::size() const {
  int size = static_cast<int>(params_[0]->size());
  for (const std::shared_ptr<ModelParam> param : params_) {
    ASSERT(size == static_cast<int>(param->size()), "size error");
  }
  return size;
}

std::shared_ptr<ModelParam> ModelParams::select(
    const std::string param_name) {
  for (std::shared_ptr<ModelParam> param : params_) {
    if (param->name() == param_name) {
      return param;
    }
  }
  ASSERT(0, "unrecognized name(" << param_name << ")");
  return NULL;
}

void ModelParams::set(const char* name,
                      const int site_type,
                      const double value) {
  select(name)->set(site_type, value);
  mix();
}

void ModelParams::set(const char* name,
    const int site_type1,
    const int site_type2,
    const double value) {
  select(name)->set_mixed(site_type1, site_type2, value);
}

}  // namespace feasst
