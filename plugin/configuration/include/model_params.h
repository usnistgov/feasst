
#ifndef FEASST_CONFIGURATION_MODEL_PARAMS_H_
#define FEASST_CONFIGURATION_MODEL_PARAMS_H_

#include <map>
#include <memory>
#include <vector>
#include <string>
#include "configuration/include/model_param.h"
#include "configuration/include/properties.h"

namespace feasst {

class Particle;
class PhysicalConstants;

typedef std::map<std::string, std::string> argtype;

/**
 The epsilon parameter is named "epsilon" in the particle file Site Types.
 The epsilon parameter has the default combining rule:
 \f$ \epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j} \f$
 */
class Epsilon : public ModelParam {
 public:
  Epsilon() : ModelParam() { class_name_ = "epsilon"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Epsilon>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<Epsilon>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Epsilon(std::istream& istr);
  virtual ~Epsilon() {}

 private:
  double mix_(const double value1, const double value2) override;
};

/**
  The sigma parameter is named "sigma" in the particle file Site Types.
  The sigma parameter has the default combining rule:

  \f$ \sigma_{ij} = \left\{
    \begin{array}{lr}
      0 & : \sigma_i \sigma_j = 0 \\
      0.5(\sigma_i + \sigma_j) & : \sigma_i\sigma_j \neq 0
    \end{array}
  \right. \f$
 */
class Sigma : public ModelParam {
 public:
  Sigma() : ModelParam() { class_name_ = "sigma"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Sigma>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<Sigma>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Sigma(std::istream& istr);
  virtual ~Sigma() {}
};

/**
  The cut off parameter is named "cutoff" in the particle file Site Types.
  The cut off parameter has the default combining rule:

  \f$ r_{c,ij} = \left\{
    \begin{array}{lr}
      0 & : r_{c,i} r_{c,j} = 0 \\
      0.5(r_{c,i} + r_{c,j}) & : r_{c,i}r_{c,j} \neq 0
    \end{array}
  \right. \f$
 */
class CutOff : public ModelParam {
 public:
  CutOff() : ModelParam() { class_name_ = "cutoff"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<CutOff>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<CutOff>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CutOff(std::istream& istr);
  virtual ~CutOff() {}
};

/**
 The charge parameter is named "charge" in the particle file Site Types.
 The charge parameter, q, has the default combining rule:
 \f$ q_{ij} = q_i q_j \f$
 */
class Charge : public ModelParam {
 public:
  Charge() : ModelParam() { class_name_ = "charge"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Charge>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<Charge>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Charge(std::istream& istr);
  virtual ~Charge() {}

 private:
  double mix_(const double value1, const double value2) override {
    return value1*value2; }
};

/**
  Container for all model parameters.
 */
class ModelParams : public PropertiedEntity {
 public:
  ModelParams();

  /// Add all site types in particle.
  void add(const Particle& particle);

  /// Compute the mixed values.
  void mix();

  /// Return the number of values.
  int size() const;

  /// Modify model parameter of a given site type and name to value.
  void set(const std::string& name, const int site_type, const double value);

  /// Modify the mixed model parameter of a pair of given site types and name
  /// to value.
  void set(const std::string& name, const int site_type1, const int site_type2,
    const double value);

  /// Set mixed model parameters by file.
  void set(const std::string& name, const std::string& filename,
    /// Optionally provide a list of names of site types
    std::vector<std::string> * site_type_names = NULL);

  /// Set mixed model parameters by file.
  void set(const std::string& filename,
    /// Optionally provide a list of names of site types
    std::vector<std::string> * site_type_names = NULL);

  /// Add a custom model parameter
  void add(std::shared_ptr<ModelParam> param) {
    params_.push_back(param); }

  /// Return the model parameter with the corresponding index.
  const ModelParam& select(const int index) const;

  /// Return the index of the model parameter with the corresponding name.
  /// Return -1 if name is not found.
  int index(const std::string& name) const;

  /// Return the model parameter with the corresponding name.
  const ModelParam& select(const std::string& name) const;

  /// Set the minimum cutoff to sigma.
  /// This is used for HardSphere potentials that don't assign cutoff.
  void set_cutoff_min_to_sigma();

  /// Set the physical constants.
  void set_physical_constants(std::shared_ptr<PhysicalConstants> constants);

  /// Set the physical constants to their default values CODATA2018.
  void set_physical_constants();

  /// Return the physical constants.
  const PhysicalConstants& physical_constants() const;

  /// Return the physical constants.
  const PhysicalConstants& constants() const { return physical_constants(); }

  /// Check
  void check() const;

  /// Return as a human readable string.
  std::string str(
     /// If provided, write the names of the site types instead of the indices.
     std::vector<std::string> * site_type_names = NULL) const;

  /// Return a deep copy of self.
  ModelParams deep_copy() const;

  void serialize(std::ostream& ostr) const;
  explicit ModelParams(std::istream& istr);

 private:
  std::vector<std::shared_ptr<ModelParam> > params_;
  std::shared_ptr<PhysicalConstants> physical_constants_;

  /// Add built-in types to params
  void add_();

  void factory_(const std::string& name);
  std::shared_ptr<ModelParam> select_(const std::string name);
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_MODEL_PARAMS_H_
