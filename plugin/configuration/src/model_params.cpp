
#include <algorithm>
#include "configuration/include/model_params.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"

namespace feasst {

void ModelParam::serialize(std::ostream& ostr) const {
  feasst_serialize(name_, ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize(values_, ostr);
  feasst_serialize(mixed_values_, ostr);
  feasst_serialize(is_mixed_override_, ostr);
}

ModelParam::ModelParam(std::istream& istr) {
  feasst_deserialize(&name_, istr);
  feasst_deserialize_version(istr);
  feasst_deserialize(&values_, istr);
  feasst_deserialize(&mixed_values_, istr);
  feasst_deserialize(&is_mixed_override_, istr);
  if (values_.size() > 0) {
    max_value_ = *std::max_element(values_.begin(), values_.end());
    max_mixed_value_ = maximum(mixed_values_);
  }
}

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
  override_resize_();
  for (int index1 = 0; index1 < size(); ++index1) {
    const double value1 = values_[index1];
    for (int index2 = 0; index2 < size(); ++index2) {
      if (!is_mixed_override_[index1][index2]) {
        const double value2 = values_[index2];
        mixed_values_[index1][index2] = mix_(value1, value2);
      }
    }
  }
  max_mixed_value_ = maximum(mixed_values_);
}

void ModelParam::set_mixed(const int site_type1,
    const int site_type2,
    const double value) {
  mixed_values_[site_type1][site_type2] = value;
  mixed_values_[site_type2][site_type1] = value;
  override_resize_();
  is_mixed_override_[site_type1][site_type2] = true;
  is_mixed_override_[site_type2][site_type1] = true;
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

void ModelParam::set_param(const ModelParams& existing) {
  // initialize the number of types with zero value
  for (int type = 0; type < existing.size(); ++type) {
    add(0.);
  }
  mix();

  // set the values
  for (int type = 0; type < existing.size(); ++type) {
    set(type, compute(type, existing));
  }

  // set the values of the mixed types
  for (int type1 = 0; type1 < existing.size(); ++type1) {
    for (int type2 = 0; type2 < existing.size(); ++type2) {
      set_mixed(type1, type2, compute(type1, type2, existing));
    }
  }
}

void ModelParam::override_resize_() {
  // assumes square size
  const int previous_size = static_cast<int>(is_mixed_override_.size());
  resize(size(), size(), &is_mixed_override_);
  for (int index1 = previous_size; index1 < size(); ++index1) {
    for (int index2 = previous_size; index2 < size(); ++index2) {
      is_mixed_override_[index1][index2] = false;
    }
  }
}

void ModelParams::add_() {
  add(epsilon_);
  add(sigma_);
  add(cutoff_);
  add(charge_);
}

ModelParams::ModelParams() : PropertiedEntity() {
  epsilon_ = std::make_shared<Epsilon>();
  sigma_ = std::make_shared<Sigma>();
  cutoff_ = std::make_shared<CutOff>();
  charge_ = std::make_shared<Charge>();
  add_();
  set_physical_constants();
}

ModelParams::ModelParams(const ModelParams& params) : PropertiedEntity() {
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
  set_physical_constants(params.physical_constants_);
}

void ModelParams::add(const Site site) {
  std::vector<std::string> names = site.properties().names();
  for (const std::string name : names) {
    std::shared_ptr<ModelParam> param = select_(name);
    if (param) {
      param->add(site);
//    } else {
//      auto param = std::make_shared<ModelParam>();
//      param->set_name(name);
//      param->add(site);
//      add(param);
    }
  }
}

void ModelParams::add(const Particle particle) {
  for (std::shared_ptr<ModelParam> param : params_) {
    param->add(particle);
//for (const Site& site : particle.sites()) {
//    add(site);
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

const std::shared_ptr<ModelParam> ModelParams::select(
    const std::string param_name) const {
  for (const std::shared_ptr<ModelParam> param : params_) {
    if (param->name() == param_name) {
      return param;
    }
  }
  ASSERT(0, "unrecognized name(" << param_name << ")");
  return NULL;
}

std::shared_ptr<ModelParam> ModelParams::select_(
    const std::string param_name) {
  for (std::shared_ptr<ModelParam> param : params_) {
    if (param->name() == param_name) {
      return param;
    }
  }
  ASSERT(0, "unrecognized name(" << param_name << ")");
  return NULL;
}

void ModelParams::set(const std::string name,
                      const int site_type,
                      const double value) {
  select_(name)->set(site_type, value);
  mix();
}

void ModelParams::set(const std::string name,
    const int site_type1,
    const int site_type2,
    const double value) {
  select_(name)->set_mixed(site_type1, site_type2, value);
}

void ModelParams::set_physical_constants(
    std::shared_ptr<PhysicalConstants> constants) {
  physical_constants_ = constants;
}

void ModelParams::check() const {
  const int size = params_.front()->size();
  for (std::shared_ptr<ModelParam> parm : params_) {
    ASSERT(size == parm->size(), "size mismatch");
  }
  properties().check();
}

// HWH warning: magic number "4"
void ModelParams::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  ostr << MAX_PRECISION;
  feasst_serialize_version(1, ostr);
  epsilon_->serialize(ostr);
  sigma_->serialize(ostr);
  cutoff_->serialize(ostr);
  charge_->serialize(ostr);
  ostr << params_.size() << " ";
  for (int index = 4; index < static_cast<int>(params_.size()); ++index) {
    params_[index]->serialize(ostr);
  }
  feasst_serialize_fstdr(physical_constants_, ostr);
}

// HWH warning: magic number "4"
ModelParams::ModelParams(std::istream& istr)
  : PropertiedEntity(istr) {
  feasst_deserialize_version(istr);
  epsilon_ = std::make_shared<Epsilon>(istr);
  sigma_ = std::make_shared<Sigma>(istr);
  cutoff_ = std::make_shared<CutOff>(istr);
  charge_ = std::make_shared<Charge>(istr);
  add_();
  int num;
  istr >> num;
  params_.resize(num);
  for (int index = 4; index < num; ++index) {
    params_[index] = std::make_shared<ModelParam>(istr);
  }
  // feasst_deserialize_fstdr(&physical_constants_, istr);
  // HWH for unknown reasons, this function template does not work.
  { int existing;
    istr >> existing;
    if (existing != 0) {
      physical_constants_ = physical_constants_->deserialize(istr);
    }
  }
}

}  // namespace feasst
