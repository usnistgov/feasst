
#ifndef FEASST_CORE_SYSTEM_H_
#define FEASST_CORE_SYSTEM_H_

#include <vector>
#include <memory>
#include <string>
#include "core/include/debug.h"
#include "core/include/configuration.h"
#include "core/include/visit_model.h"
#include "core/include/model_empty.h"
#include "core/include/constants.h"

namespace feasst {

/**
  A potential consists of a visitor, a group_index and a model.
  One issue is these potentials may need their own sets of model parameters, but these are contained
  in unique_types of configurations.
  The visitor usually takes the Model params from the configuration.
  Instead, the visitor can be told to use a different set of model parameter
 */
class Potential {
 public:
  Potential() { model_ = std::make_shared<ModelEmpty>(); }
  int group_index() const { return group_index_; }
  void set_group_index(const int group_index) {
    group_index_ = group_index;
  }
  const std::shared_ptr<VisitModel> visit_model() const {
    return visit_model_;
  }
  void set_visit_model(std::shared_ptr<VisitModel> visit_model) {
    visit_model_ = visit_model;
  }
  const std::shared_ptr<Model> model() const {
    return model_;
  }
  void set_model(std::shared_ptr<Model> model) {
    model_ = model;
  }

  double energy(Configuration * config) {
    stored_energy_ = model_->compute(group_index_, config, visit_model_.get());
    return stored_energy_;
  }

  double energy(const Select& select, Configuration * config) {
    stored_energy_ = model_->compute(select, group_index_, config, visit_model_.get());
    return stored_energy_;
  }

  double stored_energy() const { return stored_energy_; }

  void revert() { visit_model_->revert(); }

 private:
  int group_index_ = 0;
  std::shared_ptr<VisitModel> visit_model_;
  std::shared_ptr<Model> model_;
  double stored_energy_;
};

class Potentials {
 public:
  void add_potential(const Potential potential) {
    potentials_.push_back(potential); }

  double energy(Configuration * config) {
    double en = 0;
    int index = 0;
    while ((index < static_cast<int>(potentials_.size())) &&
           (en < NEAR_INFINITY)) {
      en += potentials_[index].energy(config);
      ++index;
    }
    //INFO("en " << en);
    return en;
  }

  double energy(const Select& select, Configuration * config) {
    double en = 0;
    int index = 0;
    while ((index < static_cast<int>(potentials_.size())) &&
           (en < NEAR_INFINITY)) {
      en += potentials_[index].energy(select, config);
      ++index;
    }
    //INFO("en " << en);
    return en;
  }

  void revert() {
    for (Potential& potential : potentials_) {
      potential.revert();
    }
  }

  const std::vector<Potential>& potentials() const { return potentials_; }

  double stored_energy() const {
    double en = 0.;
    for (const Potential& potential : potentials_) {
      en += potential.stored_energy();
    }
    return en;
  }

 private:
  std::vector<Potential> potentials_;
};

/**
  Systems may have multiple configurations but their typing and grouping should be the same.
  HWH refactor how the configurations are set up (e.g., no add_configuration).
  This way we can enforce typing.
  Allow duplication of configuration.
  Or maybe this should be done in the configuration class itself?
 */
class System {
 public:
  /// Set the configuration.
  void add_configuration(const Configuration& configuration) { configurations_.push_back(configuration); }

  /// Return the configuration
  const Configuration& configuration(const int index = 0) const { return configurations_[index]; }
  Configuration* get_configuration(const int index = 0) { return &configurations_[index]; }

  int dimension() const { return configurations_.front().dimension(); }

  double energy() {
    return full_.energy(&configurations_.front());
  }
  double reference_energy(const int index = 0) {
    return references_[index].energy(&configurations_.front());
  }

  double energy(const Select& select) {
    return full_.energy(select, &configurations_.front());
  }

  void set_full(const Potentials& full) { full_ = full; }
  void add_reference(const Potentials& ref) { references_.push_back(ref); }
  const Potentials& full() const { return full_; }

  void revert() { full_.revert(); }

  /// Return the header of the status of the system for periodic output.
  std::string status_header() const {
    return std::string("energy ");
  }

  /// Return the status of the system for periodic output.
  std::string status() const {
    std::stringstream ss;
    ss << full().stored_energy() << " ";
    return ss.str();
  }

 private:
  std::vector<Configuration> configurations_;

  /**
    The first potential is the full system without optimizations
    The second potential is the full system with optimizations
    The remaining potentials are used for cheap energy calculations in configurational bias
    They can also be used for reference calculations (e.g., hard sphere mayer sampling)
  */
  Potentials full_, optimized_;
  std::vector<Potentials> references_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_H_
