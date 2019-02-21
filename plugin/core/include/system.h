
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

  int dimension() const { return configurations_.front().dimension(); }

  double unoptimized_energy() {
    return unoptimized_.energy(&configurations_.front());
  }
  double energy() {
    if (is_optimized_) {
      return unoptimized_.energy(&configurations_.front());
    }
    return unoptimized_.energy(&configurations_.front());
  }
  double reference_energy(const int index = 0) {
    return references_[index].energy(&configurations_.front());
  }

  double energy(const Select& select) {
    if (is_optimized_) {
      return optimized_.energy(select, &configurations_.front());
    }
    return unoptimized_.energy(select, &configurations_.front());
  }

  void add_to_unoptimized(const Potential& potential) { unoptimized_.add_potential(potential); }
  void add_to_optimized(const Potential& potential) {
    is_optimized_ = true;
    optimized_.add_potential(potential); }
  void add_to_reference(const Potential& ref, const int index = 0) {
    if (index == 0 && references_.size() == 0) {
      references_.push_back(PotentialFactory());
    }
    references_[index].add_potential(ref);
  }
  const PotentialFactory& unoptimized() const { return unoptimized_; }
  const PotentialFactory& optimized() const { return optimized_; }

  void revert() { unoptimized_.revert(); }

  /// Return the header of the status for periodic output.
  std::string status_header() const {
    return std::string("");
  }

  /// Return the status for periodic output.
  std::string status() const {
    return std::string("");
  }

  // HWH: depreciate, but used by ewald
  void set_unoptimized(const PotentialFactory& unoptimized) { unoptimized_ = unoptimized; }

 private:
  std::vector<Configuration> configurations_;
  PotentialFactory unoptimized_;
  PotentialFactory optimized_;
  bool is_optimized_ = false;
  std::vector<PotentialFactory> references_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_H_
