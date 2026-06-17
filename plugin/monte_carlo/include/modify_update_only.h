
#ifndef FEASST_MONTE_CARLO_MODIFY_UPDATE_ONLY_H_
#define FEASST_MONTE_CARLO_MODIFY_UPDATE_ONLY_H_

#include <memory>
#include <vector>
#include <string>
#include <map>
#include "monte_carlo/include/modify.h"

namespace feasst {

class Random;
class MonteCarlo;
class TrialFactory;

/**
  This Modify does not perform writes.
  This class is for developers.
 */
class ModifyUpdateOnly : public Modify {
 public:
  explicit ModifyUpdateOnly(argtype * args);

  void set_trials_per_write(const int trials) override;

  void set_trials_per(const int trials) { set_trials_per_update(trials); }

  explicit ModifyUpdateOnly(std::istream& istr) : Modify(istr) {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_MODIFY_UPDATE_ONLY_H_
