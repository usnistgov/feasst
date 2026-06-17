
#include <cmath>
#include <algorithm>
#include "utils/include/file.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/utils.h"  // resize
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/particle.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/model_param.h"
#include "configuration/include/model_params.h"

namespace feasst {

FEASST_MAPPER(ModelParam,);

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
  set_max_and_mixed();
}

void ModelParam::set_max_and_mixed() {
  if (values_.size() > 0) {
    max_value_ = *std::max_element(values_.begin(), values_.end());
    if (mixed_values_.size() > 0) {
      max_mixed_value_ = maximum(mixed_values_);
    }
  }
}

std::map<std::string, std::shared_ptr<ModelParam> >&
    ModelParam::deserialize_map() {
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
  set_max_and_mixed();
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

std::string ModelParam::str(std::vector<std::string> * site_type_names) const {
  std::stringstream ss;
  for (int type1 = 0; type1 < size(); ++type1) {
    ss << class_name() << " ";
    if (site_type_names) {
      ss << (*site_type_names)[type1];
    } else {
      ss << type1;
    }
    ss << " " << MAX_PRECISION << value(type1) << std::endl;
  }
  for (int type1 = 0; type1 < size(); ++type1) {
    for (int type2 = type1; type2 < size(); ++type2) {
      ss << class_name() << " ";
      if (site_type_names) {
        ss << (*site_type_names)[type1] << " " << (*site_type_names)[type2];
      } else {
        ss << type1 << " " << type2;
      }
      ss << " " << MAX_PRECISION
         << mixed_value(type1, type2) << std::endl;
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

std::shared_ptr<ModelParam> ModelParam::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<ModelParam> ModelParam::factory(const std::string name,
    argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

}  // namespace feasst
