
#ifndef FEASST_MONTE_CARLO_RUN_H_
#define FEASST_MONTE_CARLO_RUN_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Perform a number of trials.
 */
class Run : public Action {
 public:
  //@{
  /** @name Arguments
    args:
    - num_trials: run this many trials (default: -1. e.g., None)
    - until_num_particles: run until this many particles (default: -1. e.g., None)
    - configuration_index: configuration for until_num_particles (default: 0).
    - particle_type: type of particle to count. If -1, all particles (default: -1).
    - for_hours: run for this many CPU hours (default: -1 e.g., None).
    - until_criteria_complete: run until Criteria is complete (default: false)

    Arguments are completed in the order listed.
   */
  explicit Run(argtype args = argtype());
  explicit Run(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Run>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Run>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Run(std::istream& istr);
  virtual ~Run() {}

  //@}
 private:
  int num_trials_;
  int until_num_particles_;
  int configuration_index_;
  int particle_type_;
  double for_hours_;
  bool until_criteria_complete_;
};

inline std::shared_ptr<Run> MakeRun(argtype args = argtype()) {
  return std::make_shared<Run>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_RUN_H_
