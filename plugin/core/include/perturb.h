
#ifndef FEASST_CORE_PERTURB_H_
#define FEASST_CORE_PERTURB_H_

#include "core/include/system.h"
#include "core/include/select_list.h"

namespace feasst {

/**
  Perturb the system (e.g., displace, add or delete particles).
  Importantly, these moves are reversible upon calling the revert function.
  Different levels of optimization govern how the revert takes place.
 */
class Perturb {
 public:
  /// Initialize some variables before each attempt.
  virtual void before_attempt() { revert_possible_ = false; }

  /// Revert the perturbation.
  virtual void revert() { *system_ = system_old_; }

  /// Return whether it is possible to revert.
  bool revert_possible() const { return revert_possible_; }

  /// Set whether it is possible to revert.
  void set_revert_possible(const bool revert = true) {
    revert_possible_ = revert;
  }

  /// Return the selection
  virtual const Select& selection() const = 0;

  virtual ~Perturb() {}

 protected:
  /* The following protected functions are only to be used by developers */

  // Before each perturbation, store the old system.
  virtual void store_old(System * system);

  System* system() { return system_; }

  int optimization_ = 1;

 private:
  System * system_;
  System system_old_;
  bool revert_possible_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_H_
