
#ifndef FEASST_CONFIGURATION_MODEL_PARAM_H_
#define FEASST_CONFIGURATION_MODEL_PARAM_H_

#include <map>
#include <memory>
#include <vector>
#include <string>
#include "configuration/include/properties.h"

namespace feasst {

class ModelParams;
class Particle;
class PhysicalConstants;
class Site;

typedef std::map<std::string, std::string> argtype;

/**
  Model parameters depend upon site types, such as epsilon, sigma, etc.
  These parameters may also be mixed between two different site types.
  Each model parameter has an assumed combining rule.
 */
class ModelParam {
 public:
  ModelParam() {}

  /// Add a new site type.
  void add(const double value);

  /// Set the value of the site type.
  void set(const int site_type, const double value) {
    values_[site_type] = value; }

  /// Add a new site type. If the site does not contain the model parameter,
  /// then add the default value.
  void add(const Site& site, const double default_value = 0);

  /// Add all site types in particle.
  void add(const Particle& particle);

  /// Compute the mixed values.
  void mix();

  /// Set the mixed value of the site types.
  void set_mixed(const int site_type1,
    const int site_type2,
    const double value);

  /// Compute and store the maximum mixed value.
  void set_max_and_mixed();

  /// Return the value.
  double value(const int type) const;

  /// Return the values.
  const std::vector<double>& values() const { return values_; }

  /// Return the mixed value.
  double mixed_value(const int type1, const int type2) const {
    return mixed_values_[type1][type2];
  }

  /// Return the mixed values.
  const std::vector<std::vector<double> >& mixed_values() const {
    return mixed_values_; }

  /// Return the number of values.
  int size() const { return static_cast<int>(values_.size()); }

  /// Return the maximum.
  double max() const { return max_value_; }

  /// Return the mixed maximum.
  double mixed_max() const;

  virtual double compute(const int type1, const int type2,
    const ModelParams& model_params) { return 0.; }

  virtual double compute(const int type1, const ModelParams& model_params) {
    return 0.; }

  /// Define new parameters modeled after the existing ones
  virtual void set_param(const ModelParams& existing);

  /// Return as a human readable string.
  std::string str(
     /// If provided, write the names of the site types instead of the indices.
     std::vector<std::string> * site_type_names = NULL) const;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<ModelParam> create(std::istream& istr) const;
  virtual std::shared_ptr<ModelParam> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<ModelParam> >& deserialize_map();
  std::shared_ptr<ModelParam> deserialize(std::istream& istr);
  std::shared_ptr<ModelParam> factory(const std::string name, argtype * args);
  explicit ModelParam(std::istream& istr);
  virtual ~ModelParam() {}

 protected:
  std::string class_name_ = "ModelParam";
  void serialize_model_param_(std::ostream& ostr) const;

 private:
  std::vector<double> values_;
  std::vector<std::vector<double> > mixed_values_;
  double max_value_;
  double max_mixed_value_;
  std::vector<std::vector<bool> > is_mixed_override_;

  /// Define combining rules in the derived class.
  /// The default is a simple average, unless one of the values is zero.
  virtual double mix_(const double value1, const double value2);

  // resize is_mixed_override_ by setting to false by default but not
  // overwriting any that were previously set to true.
  void override_resize_();
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_MODEL_PARAM_H_
