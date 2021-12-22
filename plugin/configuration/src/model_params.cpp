
#include <cmath>
#include <algorithm>
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "utils/include/utils.h"  // resize
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"

namespace feasst {

class MapModelParam {
 public:
  MapModelParam() {
    auto obj = std::make_shared<ModelParam>();
    obj->deserialize_map()["ModelParam"] = obj;
  }
};

static MapModelParam mapper_ = MapModelParam();

void ModelParam::serialize_model_param_(std::ostream& ostr) const {
  feasst_serialize_version(4795, ostr);
  feasst_serialize(values_, ostr);
  feasst_serialize(mixed_values_, ostr);
  feasst_serialize(is_mixed_override_, ostr);
}

void ModelParam::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  serialize_model_param_(ostr);
}

ModelParam::ModelParam(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4795, "unrecognized version: " << version);
  feasst_deserialize(&values_, istr);
  feasst_deserialize(&mixed_values_, istr);
  feasst_deserialize(&is_mixed_override_, istr);
  if (values_.size() > 0) {
    max_value_ = *std::max_element(values_.begin(), values_.end());
    if (mixed_values_.size() > 0) {
      max_mixed_value_ = maximum(mixed_values_);
    }
  }
}

std::map<std::string, std::shared_ptr<ModelParam> >& ModelParam::deserialize_map() {
  static std::map<std::string, std::shared_ptr<ModelParam> >* ans =
     new std::map<std::string, std::shared_ptr<ModelParam> >();
  return *ans;
}

std::shared_ptr<ModelParam> ModelParam::create(std::istream& istr) const {
  return std::make_shared<ModelParam>(istr);
}

std::shared_ptr<ModelParam> ModelParam::create(argtype * args) const {
  FATAL("not implemented");
}

void ModelParam::add(const double value) {
  values_.push_back(value);
  max_value_ = *std::max_element(values_.begin(), values_.end());
}

void ModelParam::add(const Site& site, const double default_value) {
  DEBUG(class_name() << " adding site");
  double value;
  if (site.properties().value(class_name(), &value)) {
    add(value);
  } else {
    add(default_value);
  }
}

void ModelParam::add(const Particle& particle) {
  DEBUG("adding");
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
  for (int type = size(); type < existing.size(); ++type) {
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

std::string ModelParam::str() const {
  std::stringstream ss;
  ss << class_name() << std::endl;
  for (int type1 = 0; type1 < size(); ++type1) {
    for (int type2 = 0; type2 < size(); ++type2) {
      ss << "type" << type1 << "," << value(type1) << ",";
      ss << "type" << type2 << "," << value(type2) << ",";
      if (type1 < static_cast<int>(mixed_values().size())) {
        if (type2 < static_cast<int>(mixed_values()[type1].size())) {
          ss << "type" << type1 << "-" << type2 << "," << mixed_value(type1, type2);
        }
      }
      ss << std::endl;
    }
  }
  return ss.str();
}

double ModelParam::mix_(const double value1, const double value2) {
  if (std::abs(value1) < NEAR_ZERO ||
      std::abs(value2) < NEAR_ZERO) {
    return 0.;
  }
  return 0.5*(value1 + value2);
}

ModelParams::ModelParams() : PropertiedEntity() {
//  add(std::make_shared<Epsilon>());
//  add(std::make_shared<Sigma>());
//  add(std::make_shared<CutOff>());
//  add(std::make_shared<Charge>());
  set_physical_constants();
}

//void ModelParams::add(const Site site) {
//  std::vector<std::string> names = site.properties().names();
//  DEBUG("adding: " << feasst_str(names));
//  for (const std::string name : names) {
//    DEBUG("name: " << name);
//    std::shared_ptr<ModelParam> param = select_(name);
//    if (!param) {
//      // search for param in deserialize_map and add it
//      param = ModelParam().factory(name, NULL);
//    }
//    if (param) {
//      param->add(site);
////    } else {
////      auto param = std::make_shared<ModelParam>();
////      param->set_name(name);
////      param->add(site);
////      add(param);
//    }
//  }
//}

void ModelParams::add(const Particle& particle) {
  DEBUG("adding");
  // new params added based on properties of unique site types
  for (const Site& site : particle.sites()) {
    for (const std::string name : site.properties().names()) {
      if (ModelParam().deserialize_map().count(name) != 0) {
        bool found = false;
        for (std::shared_ptr<ModelParam> param : params_) {
          if (name == param->class_name()) {
            found = true;
          }
        }
        if (!found) {
          DEBUG("automatically adding " << name);
          add(deep_copy_derived(ModelParam().deserialize_map()[name]));
        }
      }
    }
  }

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
  if (params_.size() == 0) {
    return 0;
  }
  int size = static_cast<int>(params_[0]->size());
  for (const std::shared_ptr<ModelParam> param : params_) {
    ASSERT(size == static_cast<int>(param->size()), "size error");
  }
  return size;
}

int ModelParams::index(const std::string param_name) const {
  for (int ind = 0; ind < static_cast<int>(params_.size()); ++ind) {
    const std::shared_ptr<ModelParam> param = params_[ind];
    ASSERT(param, "error");
    if (param->class_name() == param_name) {
      return ind;
    }
  }
  return -1;
  FATAL("unrecognized name(" << param_name << ")");
}

const ModelParam& ModelParams::select(const int index) const {
  ASSERT(index >= 0, "index: " << index);
  ASSERT(index < static_cast<int>(params_.size()),
    "index: " << index << " >= size: " << params_.size());
  return const_cast<const ModelParam&>(*params_[index]);
}

const ModelParam& ModelParams::select(const std::string name) const {
  const int indx = index(name);
  ASSERT(indx != -1, "name: " << name << " not found");
  return const_cast<const ModelParam&>(*params_[indx]);
}

std::shared_ptr<ModelParam> ModelParams::select_(
    const std::string param_name) {
  for (std::shared_ptr<ModelParam> param : params_) {
    if (param->class_name() == param_name) {
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

void ModelParams::serialize(std::ostream& ostr) const {
  PropertiedEntity::serialize(ostr);
  ostr << MAX_PRECISION;
  feasst_serialize_version(938, ostr);
  feasst_serialize(params_, ostr);
//  ostr << params_.size() << " ";
//  for (int index = 0; index < static_cast<int>(params_.size()); ++index) {
//    params_[index]->serialize(ostr);
//  }
  feasst_serialize_fstdr(physical_constants_, ostr);
}

ModelParams::ModelParams(std::istream& istr)
  : PropertiedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 938, "unrecognized version: " << version);
//  int num;
//  istr >> num;
//  params_.resize(num);
//  for (int index = 0; index < num; ++index) {
//    params_[index] = std::make_shared<ModelParam>(istr);
//  }
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize(&params_, istr);
  { int dim1;
    istr >> dim1;
    params_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstobj((params_)[index], istr);
      int existing;
      istr >> existing;
      DEBUG("existing: " << existing);
      if (existing != 0) {
        params_[index] = params_[index]->deserialize(istr);
      }
    }
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

double Epsilon::mix_(const double value1, const double value2) {
  return std::sqrt(value1*value2);
}

void ModelParams::set_cutoff_min_to_sigma() {
  DEBUG("here goes");
  auto sigma = select_("sigma");
  auto cutoff = select_("cutoff");
  DEBUG("here goes");
  for (int itype = 0; itype < sigma->size(); ++itype) {
    for (int jtype = 0; jtype < sigma->size(); ++jtype) {
      const double sig = sigma->mixed_value(itype, jtype);
      if (cutoff->mixed_value(itype, jtype) < sig) {
        cutoff->set_mixed(itype, jtype, sig);
      }
    }
  }
}

std::string ModelParams::str() const {
  std::stringstream ss;
  for (const std::shared_ptr<ModelParam> param : params_) {
    ss << param->str();
  }
  return ss.str();
}

//ModelParams::ModelParams(const ModelParams& params) : PropertiedEntity() {
//  FATAL("ModelParams copy constructor seems to have issues with shared_ptr");
//  for (int param = 0; param < static_cast<int>(params.params_.size());
//       ++param) {
//    add(params.params_[param]);
//  }
//  set_properties(params.properties());
//  set_physical_constants(params.physical_constants_);
//}

ModelParams ModelParams::deep_copy() const {
  return feasst::deep_copy(*this);
}

std::shared_ptr<ModelParam> ModelParam::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<ModelParam> ModelParam::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

class MapEpsilon {
 public:
  MapEpsilon() {
    auto obj = std::make_shared<Epsilon>();
    obj->deserialize_map()["epsilon"] = obj;
  }
};

static MapEpsilon mapper_epsilon_ = MapEpsilon();

void Epsilon::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2048, ostr);
}

Epsilon::Epsilon(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2048, "mismatch version: " << version);
}

class MapSigma {
 public:
  MapSigma() {
    auto obj = std::make_shared<Sigma>();
    obj->deserialize_map()["sigma"] = obj;
  }
};

static MapSigma mapper_sigma_ = MapSigma();

void Sigma::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(6095, ostr);
}

Sigma::Sigma(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6095, "mismatch version: " << version);
}

class MapCutOff {
 public:
  MapCutOff() {
    auto obj = std::make_shared<CutOff>();
    obj->deserialize_map()["cutoff"] = obj;
  }
};

static MapCutOff mapper_cutoff_ = MapCutOff();

void CutOff::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2948, ostr);
}

CutOff::CutOff(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2948, "mismatch version: " << version);
}

class MapCharge {
 public:
  MapCharge() {
    auto obj = std::make_shared<Charge>();
    obj->deserialize_map()["charge"] = obj;
  }
};

static MapCharge mapper_charge_ = MapCharge();

void Charge::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(1094, ostr);
}

Charge::Charge(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1094, "mismatch version: " << version);
}

}  // namespace feasst
