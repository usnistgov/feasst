
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

ModelParams::ModelParams() {
  epsilon_ = std::make_shared<Epsilon>();
  sigma_ = std::make_shared<Sigma>();
  cutoff_ = std::make_shared<CutOff>();
  charge_ = std::make_shared<Charge>();
  params_.push_back(epsilon_);
  params_.push_back(sigma_);
  params_.push_back(cutoff_);
  params_.push_back(charge_);
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

std::shared_ptr<ModelParam> ModelParams::selector_(
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
  selector_(name)->set(site_type, value);
  mix();
}

}  // namespace feasst
