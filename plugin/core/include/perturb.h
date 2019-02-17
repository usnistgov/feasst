
#ifndef FEASST_CORE_PERTURB_H_
#define FEASST_CORE_PERTURB_H_

#include "core/include/system.h"
#include "core/include/select_list.h"

namespace feasst {

/**
  Perturb the system (e.g., displace, add or delete particles).
  Importantly, these moves are reversible upon calling the revert function.
 */
class Perturb {
 public:
  /// Initialize some variables before each attempt.
  virtual void before_attempt() { revert_possible_ = false; }

  /// Revert the perturbation.
  virtual void revert() {
    ASSERT(!optimized_revert(), "nonoptimized revert requires system storage");
    *system_ = system_old_;
  }

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
  void store_old(System * system) {
    system_ = system;
    if (!optimized_revert()) {
      system_old_ = *system;
    }
  }

  virtual bool optimized_revert() { return false; }

  System* system() { return system_; }

 private:
  System * system_;
  System system_old_;
  bool revert_possible_;
};

class PerturbOptRevert : public Perturb {
  bool optimized_revert() override { return true; }
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_H_
