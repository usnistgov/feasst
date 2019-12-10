
#ifndef FEASST_SYSTEM_POTENTIAL_H_
#define FEASST_SYSTEM_POTENTIAL_H_

#include <memory>
#include "utils/include/cache.h"
#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"
#include "system/include/model.h"

namespace feasst {

/**
  A potential represents both the model and the method used to compute the
  model (e.g., VisitModel) for a given Configuration.
  A potential also allows customization of the model parameters templated
  but separate from those in the Configuration.
 */
class Potential {
 public:
  /**
    Construct with the default visitor and default model (emtpy).

    args:
    - group_index: set the index of the group in the configuration which
      contributes to this potential (default: 0, representing entire config).
    - cell_index: set the index of the cell, only to be used with
      VisitModelCell.
      This also overrides group_index.
   */
  explicit Potential(const argtype& args = argtype());

  /// Return the index of the group
  int group_index() const { return group_index_; }

  /// Return the index of the cell
  int cell_index() const;

  /// Constract with model and default visitor.
  explicit Potential(std::shared_ptr<Model> model,
                     const argtype& args = argtype());

  /// Return the model.
  const Model * model() const { return model_.get(); }

  /// Construct with visitor and default model.
  explicit Potential(std::shared_ptr<VisitModel> visit_model,
                     const argtype& args = argtype());

  /// Return the method used to compute.
  const VisitModel * visit_model() const { return visit_model_.get(); }

  /// Construct with model and visitor.
  Potential(std::shared_ptr<Model> model,
            std::shared_ptr<VisitModel> visit_model,
            const argtype& args = argtype());

  /// Set the model parameters. If not set, use the one from configuration.
  void set_model_params(const ModelParams& model_params) {
    model_params_override_ = true;
    model_params_ = model_params;
  }

  /// Set the model parameters to the one in the configuration.
  void set_model_params(const Configuration& config) {
    set_model_params(config.model_params()); }

  /// Modify model parameter of a given site type and name to value.
  void set_model_param(const char* name,
                       const int site_type,
                       const double value);

  /// Return the model parameters.
  const ModelParams& model_params() const;

  /// Return the model parameters.
  /// Use model parameters from configuration if they have not been overriden.
  const ModelParams& model_params(const Configuration * config) const;

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute(Configuration * config) {
    visit_model_->precompute(config);
    model_->precompute(model_params(config));
  }

  /// Compute the energy of the entire configuration.
  double energy(Configuration * config);

  /// Compute the energy of a selection of the configuration.
  double energy(const Select& select, Configuration * config);

  /// Return the last computed value of the energy.
  double stored_energy() const { return stored_energy_; }

  /// Revert any changes to the configuration due to the last energy computation
  void revert() { visit_model_->revert(); }

  /// Finalize changes to the configuration due to the last energy computation
  void finalize() { visit_model_->finalize(); }

  /// Return the cache.
  const Cache& cache() const { return cache_; }

  /// Set Cache to load.
  void load_cache(const bool load) { cache_.set_load(load); }

  /// Set Cache to unload.
  void unload_cache(const Potential& potential) {
    cache_.set_unload(potential.cache()); }

  /// Serialize.
  void serialize(std::ostream& ostr) const;

  /// Deserialize.
  explicit Potential(std::istream& istr);

 private:
  int group_index_;
  std::shared_ptr<VisitModel> visit_model_;
  std::shared_ptr<Model> model_;
  double stored_energy_ = 0.;
  bool model_params_override_ = false;
  ModelParams model_params_;
  Cache cache_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_POTENTIAL_H_
