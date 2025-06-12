
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

ModelParams::ModelParams() : PropertiedEntity() {
//  add(std::make_shared<Epsilon>());
//  add(std::make_shared<Sigma>());
//  add(std::make_shared<CutOff>());
//  add(std::make_shared<Charge>());
  set_physical_constants();
}

// void ModelParams::add(const Site site) {
//   std::vector<std::string> names = site.properties().names();
//   DEBUG("adding: " << feasst_str(names));
//   for (const std::string name : names) {
//     DEBUG("name: " << name);
//     std::shared_ptr<ModelParam> param = select_(name);
//     if (!param) {
//       // search for param in deserialize_map and add it
//       param = ModelParam().factory(name, NULL);
//     }
//     if (param) {
//       param->add(site);
// //    } else {
// //      auto param = std::make_shared<ModelParam>();
// //      param->set_name(name);
// //      param->add(site);
// //      add(param);
//     }
//   }
// }

void ModelParams::add(const Particle& particle) {
  DEBUG("adding");
  // new params added based on properties of unique site types
  for (const Site& site : particle.sites()) {
    for (const std::string& name : site.properties().names()) {
      if (ModelParam().deserialize_map().count(name) != 0) {
        bool found = false;
        for (std::shared_ptr<ModelParam> param : params_) {
          if (name == param->class_name()) {
            found = true;
          }
        }
        if (!found) {
          DEBUG("automatically adding " << name);
          auto model_param =
            deep_copy_derived(ModelParam().deserialize_map()[name]);
          // DEBUG(particle.num_sites());
          // DEBUG(size());
          for (int i = 0; i < size(); ++i) {
            model_param->add(0.);  // add the default placeholder value
          }
          add(model_param);
        }
      }
    }
  }

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
  if (params_.size() == 0) {
    return 0;
  }
  int size = static_cast<int>(params_[0]->size());
  for (const std::shared_ptr<ModelParam>& param : params_) {
    ASSERT(size == static_cast<int>(param->size()), "size error");
  }
  return size;
}

int ModelParams::index(const std::string& param_name) const {
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
  ASSERT(index >= 0, "The expected ModelParam was not found. "
    << "Check that the particle_type includes the parameters expected by "
    << "the given models. index: " << index);
  ASSERT(index < static_cast<int>(params_.size()),
    "index: " << index << " >= size: " << params_.size());
  return const_cast<const ModelParam&>(*params_[index]);
}

const ModelParam& ModelParams::select(const std::string& name) const {
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
  ASSERT(0, "Unrecognized ModelParam name(" << param_name << ")");
  return NULL;
}

void ModelParams::set(const std::string& name,
                      const int site_type,
                      const double value) {
  select_(name)->set(site_type, value);
  mix();
}

void ModelParams::set(const std::string& name,
    const int site_type1,
    const int site_type2,
    const double value) {
  select_(name)->set_mixed(site_type1, site_type2, value);
}

void ModelParams::set(const std::string& name, const std::string& filename,
    std::vector<std::string> * site_type_names) {
  WARN("Deprecated.");
  std::ifstream file(filename.c_str());
  ASSERT(file.good(), "cannot find mixing file " << filename.c_str());
  std::shared_ptr<ModelParam> param = select_(name);
  const int size = param->size();
  int site_type1, site_type2;
  std::string stname1, stname2;
  double value;
  std::string line;
  bool read = true, found;
  while (read) {
    if (site_type_names) {
      file >> stname1 >> stname2 >> value;
      found = find_in_list(stname1, *site_type_names, &site_type1);
      ASSERT(found, "Could not find site type name:" << stname1);
      found = find_in_list(stname2, *site_type_names, &site_type2);
      ASSERT(found, "Could not find site type name:" << stname2);
    } else {
      file >> site_type1 >> site_type2 >> value;
    }
    DEBUG("site_type1 " << site_type1 << " site_type2 " << site_type2 <<
      " value " << value << " size " << size);
    ASSERT(site_type1 < size && site_type2 < size,
      "given site_type1:" << site_type1 << " or site_type2: " << site_type2
      << " >= size:" << size << ". The mixing file seems to have more sites"
      << " than the current configuration.");
    param->set_mixed(site_type1, site_type2, value);
    std::getline(file, line);
    if (file.peek() == EOF) {
      read = false;
    }
  }
}

void ModelParams::set(const std::string& filename,
    std::vector<std::string> * site_type_names) {
  std::ifstream file(filename.c_str());
  ASSERT(file.good(), "cannot find mixing file " << filename.c_str());
  std::string line;
  //std::getline(file, line);
  bool last_line = false, found;
  int num = 0;
  int site_type1, site_type2;
  std::shared_ptr<ModelParam> param;
  while (!last_line) {
    if (file.eof()) last_line = true;
    std::getline(file, line);
    ++num;
    ASSERT(num < 1e9, "The file:" << filename << " has too many lines " <<
      "or is otherwise improperly formatted.");
    DEBUG("line:" << line);
    if (line.front() != '#' && !line.empty()) {
      DEBUG("parsing");
      std::vector<std::string> vals = split(line, ' ');
      param = select_(vals[0]);
      const int size = param->size();
      if (static_cast<int>(vals.size()) == 3) {
        if (site_type_names) {
          found = find_in_list(vals[1], *site_type_names, &site_type1);
          ASSERT(found, "Could not find site type name:" << vals[1]);
        } else {
          site_type1 = feasst::str_to_int(vals[1]);
        }
        ASSERT(site_type1 < size, "given site_type1:" << site_type1 <<
          " >= size:" << size << ". The file:"<<filename<<" seems to have more "
          << "sites than the current configuration:"<<param->size());
        param->set(site_type1, feasst::str_to_double(vals[2]));
      } else if (static_cast<int>(vals.size()) == 4) {
        if (site_type_names) {
          found = find_in_list(vals[1], *site_type_names, &site_type1);
          ASSERT(found, "Could not find site type name:" << vals[1]);
          found = find_in_list(vals[2], *site_type_names, &site_type2);
          ASSERT(found, "Could not find site type name:" << vals[2]);
        } else {
          site_type1 = feasst::str_to_int(vals[1]);
          site_type2 = feasst::str_to_int(vals[2]);
        }
        ASSERT(site_type1 < size && site_type2 < size,
          "given site_type1:" << site_type1 << " or site_type2: " << site_type2
          << " >= size:" << size << ". The file:"<<filename<<" seems to have "
          << "more sites than the current configuration:"<<param->size());
        param->set_mixed(site_type1, site_type2, feasst::str_to_double(vals[3]));
      } else {
        FATAL("Unrecognized ModelParams format for file:" << filename <<
          "on line:\"" << line << "\"");
      }
    }
  }
}

void ModelParams::set_physical_constants(
    std::shared_ptr<PhysicalConstants> constants) {
  physical_constants_ = constants;
}

void ModelParams::check() const {
  const int size = params_.front()->size();
  for (std::shared_ptr<ModelParam> parm : params_) {
    ASSERT(size == parm->size(), "size mismatch for " << parm->class_name() <<
      " of " << size << " vs " << parm->size());
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
  auto sigma = select_("sigma");
  auto cutoff = select_("cutoff");
  for (int itype = 0; itype < sigma->size(); ++itype) {
    const double isig = sigma->value(itype);
    DEBUG("isig " << isig);
    if (cutoff->value(itype) < isig) {
      DEBUG("setting " << itype << " sig " << isig);
      cutoff->set(itype, isig);
    }
    for (int jtype = 0; jtype < sigma->size(); ++jtype) {
      const double sig = sigma->mixed_value(itype, jtype);
      if (cutoff->mixed_value(itype, jtype) < sig) {
        DEBUG("setting " << itype << "-" << jtype << " sig " << sig);
        cutoff->set_mixed(itype, jtype, sig);
      }
    }
  }
  cutoff->set_max_and_mixed();
}

std::string ModelParams::str(std::vector<std::string> * site_type_names) const {
  std::stringstream ss;
  ss << "# This is a FEASST Configuration::model_param_file which may " <<
    "override the default Lorentz-Berthelot mixing rules." << std::endl <<
    "# Each line has a 3- or 4-column space-separated format:" << std::endl <<
    "# [parameter] [site type] [value]" << std::endl <<
    "# [parameter] [site type 1] [site type 2] [value]" << std::endl << std::endl;
  for (const std::shared_ptr<ModelParam>& param : params_) {
    ss << param->str(site_type_names);
  }
  return ss.str();
}

// ModelParams::ModelParams(const ModelParams& params) : PropertiedEntity() {
//   FATAL("ModelParams copy constructor seems to have issues with shared_ptr");
//   for (int param = 0; param < static_cast<int>(params.params_.size());
//        ++param) {
//     add(params.params_[param]);
//   }
//   set_properties(params.properties());
//   set_physical_constants(params.physical_constants_);
// }

ModelParams ModelParams::deep_copy() const {
  return feasst::deep_copy(*this);
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

FEASST_MAPPER_RENAME(Epsilon, epsilon,);

void Epsilon::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2048, ostr);
}

Epsilon::Epsilon(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2048, "mismatch version: " << version);
}

FEASST_MAPPER_RENAME(Sigma, sigma,);

void Sigma::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(6095, ostr);
}

Sigma::Sigma(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6095, "mismatch version: " << version);
}

FEASST_MAPPER_RENAME(CutOff, cutoff,);

void CutOff::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2948, ostr);
}

CutOff::CutOff(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2948, "mismatch version: " << version);
}

FEASST_MAPPER_RENAME(Charge, charge,);

void Charge::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(1094, ostr);
}

Charge::Charge(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1094, "mismatch version: " << version);
}

const PhysicalConstants& ModelParams::physical_constants() const {
  return const_cast<PhysicalConstants&>(*physical_constants_);
}

void ModelParams::set_physical_constants() {
  set_physical_constants(MakeCODATA2018());
}

}  // namespace feasst
