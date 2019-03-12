
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
#include "core/include/potential_factory.h"

namespace feasst {

/**
  Systems may have multiple configurations but their typing and grouping should be the same.
  HWH refactor how the configurations are set up (e.g., no add_configuration).
  This way we can enforce typing.
  Allow duplication of configuration.
  Or maybe this should be done in the configuration class itself?

    The first potential is the system without optimizations
    The second potential is the system with optimizations
    The remaining potentials are used for cheap energy calculations in configurational bias
    They can also be used for reference calculations (e.g., hard sphere mayer sampling)
 */
class System {
 public:
  /// Set the configuration.
  void add(const Configuration& configuration) { configurations_.push_back(configuration); }

  /// Return the configuration
  const Configuration& configuration(const int index = 0) const { return configurations_[index]; }
  Configuration* get_configuration(const int index = 0) { return &configurations_[index]; }

  int dimension(const int config_index = 0) const {
    const int dim = configurations_[0].dimension();
    for (const Configuration& config : configurations_) {
      ASSERT(dim == config.dimension(), "dimensions of configs do not match");
    }
    return dim;
  }

  double unoptimized_energy() {
    return unoptimized_.energy(&configurations_.front());
  }
  double energy() {
    // for when bonded energies, etc come into play,
    // perhaps have energy of full system use reference potentials.
    // this way CB could distinguish external and internal interactions.
    if (is_optimized_) {
      return optimized_.energy(&configurations_.front());
    }
    return unoptimized_.energy(&configurations_.front());
  }
  double reference_energy(const int index = 0) {
    return reference_(index)->energy(&configurations_.front());
  }
  double reference_energy(const Select& select, const int index = 0) {
    return reference_(index)->energy(select, &configurations_.front());
  }


  double energy(const Select& select) {
    if (is_optimized_) {
      return optimized_.energy(select, &configurations_.front());
    }
    return unoptimized_.energy(select, &configurations_.front());
  }

  void add(const Potential& potential) { add_to_unoptimized(potential); }
  void add_to_unoptimized(const Potential& potential) { unoptimized_.add_potential(potential); }
  void add_to_optimized(const Potential& potential) {
    is_optimized_ = true;
    optimized_.add_potential(potential); }
  void add_to_reference(const Potential& ref, const int index = 0) {
    if (index == 0 && references_.size() == 0) {
      references_.push_back(PotentialFactory());
    }
    reference_(index)->add_potential(ref);
  }
  const PotentialFactory& unoptimized() const { return unoptimized_; }
  const PotentialFactory& optimized() const { return optimized_; }

  void revert() { unoptimized_.revert(); }

  // HWH: depreciate, but used by ewald
  void set_unoptimized(const PotentialFactory& unoptimized) { unoptimized_ = unoptimized; }

  void precompute() {
    unoptimized_.precompute(&configurations_.front());
    if (is_optimized_) {
      optimized_.precompute(&configurations_.front());
    }
    for (PotentialFactory& ref : references_) {
      ref.precompute(&configurations_.front());
    }
  }

  const std::vector<PotentialFactory> references() const { return references_; }
  const Potential& reference(const int ref, const int potential) const {
    return references_[ref].potentials()[potential]; }

 private:
  std::vector<Configuration> configurations_;
  PotentialFactory unoptimized_;
  PotentialFactory optimized_;
  bool is_optimized_ = false;
  std::vector<PotentialFactory> references_;

  PotentialFactory * reference_(const int index) {
    ASSERT(index < static_cast<int>(references_.size()),
      "unrecognized reference");
    return &references_[index];
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_H_
