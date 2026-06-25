
#ifndef FEASST_SYSTEM_POTENTIAL_H_
#define FEASST_SYSTEM_POTENTIAL_H_

#include <map>
#include <string>
#include <memory>

namespace feasst {

class Cache;
class Configuration;
class Model;
class ModelParams;
class Select;
class VisitModel;

typedef std::map<std::string, std::string> argtype;

/**
  A Potential represents both the Model and the method used to compute the
  Model (e.g., VisitModel) for a given Configuration.
  A Potential also allows customization of the ModelParams templated
  but separate from those in the Configuration.
 */
class Potential {
 public:
  //@{
  /* Deprecated argument:
    - group_index: set the index of the Group in the Configuration which
      contributes to this Potential (default: 0, representing entire config).
      Cannot be used together with group.
    */
  /** @name Arguments
    - Model: Derived class Model name (default: ModelEmpty).
    - VisitModel: Derived class VisitModel name (default: VisitModel).
    - config: Name of Configuration in System (default: 0).
    - group: Name of Group defined within Configuration
      (default: the whole Configuration).
    - cell_index: Optionally set the index of the cell, only used with
      VisitModelCell. This also overrides group.
    - prevent_cache: Set this to true in order to prevent the use of cache
      (default: false)
    - table_size: Set size of ModelTwoBodyTable (default: 0).
      Do not use if size <= 0.
    - table_hard_sphere_threshold: If using a table above, set the
      ModelTwoBodyTable hard_sphere_threshold (default: 0.85).
    - [parameter]: Optionally, override ModelParams as described in
      Configuration arguments.
    - [parameter][type1]: Optionally, override ModelParams as described in
      Configuration arguments.
    - [parameter][type1]_[type2]: Optionally, override ModelParams as described
      in Configuration arguments.
   */
  explicit Potential(argtype args = argtype());
  explicit Potential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the index of the group
  int group_index() const { return group_index_; }

  /// Return the index of the cell
  int cell_index() const;

  /// Return the configuration name
  const std::string& config() const { return config_; }

  /// Return the configuration index
  int configuration_index() const { return configuration_index_; }

  /// Set the configuration index
  void set_configuration_index(const int configuration_index) {
    configuration_index_ = configuration_index; }

  /// Construct with model and default visitor.
  explicit Potential(std::shared_ptr<Model> model, argtype args = argtype());

  /// Return the model.
  const Model& model() const;

  /// Construct with visitor and default model.
  Potential(std::shared_ptr<VisitModel> visit_model,
            argtype args = argtype());

  /// Return the method used to compute.
  const VisitModel& visit_model() const;

  /// Construct with model and visitor.
  Potential(std::shared_ptr<Model> model,
            std::shared_ptr<VisitModel> visit_model,
            argtype args = argtype());

  /// Set the model parameters. If not set, use the one from configuration.
  void set(const ModelParams& model_params);

  bool are_model_params_overridden() const { return model_params_override_; }

  /// Set the model parameters to the one in the configuration.
  void set_model_params(const Configuration& config);

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const std::string& name,
                       const int site_type,
                       const double value);
  void set_model_param(const std::string& name,
                       const int site_type,
                       const double value,
                       const Configuration& config);

  /// Modify model parameter of given site types and name to value.
  void set_model_param(const std::string& name,
                       const int site_type0,
                       const int site_type1,
                       const double value);
  void set_model_param(const std::string& name,
                       const int site_type0,
                       const int site_type1,
                       const double value,
                       const Configuration& config);

  /// Return the model parameters.
  const ModelParams& model_params() const;

  /// Return the model parameters.
  /// Use model parameters from configuration if they have not been overridden.
  const ModelParams& model_params(const Configuration& config) const;

  /// Return the model parameters.
  /// Use model parameters from configuration if they have not been overridden.
  ModelParams * get_model_params(Configuration * config);

  /// Check that the cutoff is within the allowed range for the Domain.
  bool does_cutoff_fit_domain(const Configuration& config,
                              /// Generate a fatal error if returning false
                              const bool fatal = false) const;

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute(Configuration * config);

  /// Compute the energy of the entire configuration.
  virtual double energy(Configuration * config);

  /// Compute the energy of a selection of the configuration.
  virtual double select_energy(const Select& select, Configuration * config);

  /// Return the last computed value of the energy.
  double stored_energy() const { return stored_energy_; }

  /// Set the last computed value of the energy.
  void set_stored_energy(const double energy) { stored_energy_ = energy; }

  /// Change the volume.
  void change_volume(const double delta_volume, const int dimension,
      Configuration * config);

  /// Revert any changes to the configuration due to the last energy computation
  void revert(const Select& select);

  /// Finalize changes to the configuration due to the last energy computation
  void finalize(const Select& select, Configuration * config);

  /// Return the cache.
  const Cache& cache() const;

  /// Set Cache to load.
  void load_cache(const bool load);

  /// Set Cache to unload.
  void unload_cache(const Potential& potential);

  void synchronize_(const Potential& potential, const Select& perturbed);

  void check(const Configuration& config) const;

  void set_visit_model_(std::shared_ptr<VisitModel> visit);

  void set_model_index(const int index);

  VisitModel * get_visit_model_() { return visit_model_.get(); }

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Deserialize.
  explicit Potential(std::istream& istr);
  virtual ~Potential();

  //@}
 private:
  int group_index_;
  std::string group_;
  std::shared_ptr<VisitModel> visit_model_;
  std::shared_ptr<Model> model_;
  double stored_energy_ = 0.;
  bool model_params_override_ = false;
  std::shared_ptr<ModelParams> model_params_;
  std::shared_ptr<Cache> cache_;
  bool prevent_cache_;
  int table_size_;
  double table_hs_threshold_;
  argtype override_args_;
  int configuration_index_;
  std::string config_;
};

inline std::shared_ptr<Potential> MakePotential(argtype args = argtype()) {
  return std::make_shared<Potential>(args);
}

inline std::shared_ptr<Potential> MakePotential(
    std::shared_ptr<Model> model,
    argtype args = argtype()) {
  return std::make_shared<Potential>(model, args);
}

inline std::shared_ptr<Potential> MakePotential(
    std::shared_ptr<VisitModel> visit_model,
    argtype args = argtype()) {
  return std::make_shared<Potential>(visit_model, args);
}

inline std::shared_ptr<Potential> MakePotential(
    std::shared_ptr<Model> model,
    std::shared_ptr<VisitModel> visit_model,
    argtype args = argtype()) {
  return std::make_shared<Potential>(model, visit_model, args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_POTENTIAL_H_
