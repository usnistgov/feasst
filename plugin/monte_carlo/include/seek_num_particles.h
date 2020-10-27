
#ifndef FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_
#define FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

class Criteria;
class Constraint;
class TrialAdd;
class Trial;
class MonteCarlo;
class ProgressReport;

/**
  Initialize the number of particles in MonteCarlo.
  Use the existing trials, but reset their stats once the goal is reached.
 */
class SeekNumParticles {
 public:
  /**
    args:
    - particle_type: type of particle to add. If -1, all types (default: -1).
    - max_attempts: maximum number of trial attempts for seek (default: 1e8).
   */
  SeekNumParticles(const int num, const argtype& args = argtype());

  /// Optionally, use a temporary ThermoParams.
  SeekNumParticles with_thermo_params(const argtype& args);

  /// Optionally use a temporary Metropolis criteria with given args.
  SeekNumParticles with_metropolis();

  /// Same as above but with a constraint.
  SeekNumParticles with_metropolis(std::shared_ptr<Constraint> constraint);

  /**
    Use an extra and temporary TrialAdd between every existing Trial.
    If particle_type is not included in args, use value given by class
    argument particle_type, if not -1.
    If -1, then use 0.
   */
  SeekNumParticles with_trial_add(const argtype& args = argtype());

  /**
    Add an extra trial.
    Any extra trials, including TrialAdd above, are attempted between each trial
    already present in MonteCarlo.
   */
  SeekNumParticles add(std::shared_ptr<Trial> trial);

  /// Add a progress report.
  SeekNumParticles add(std::shared_ptr<ProgressReport> report) {
    report_ = report; return *this; }

  /// Perform the seek.
  void run(MonteCarlo * monte_carlo);

 private:
  const int num_;
  int max_attempts_;
  int particle_type_;
  std::shared_ptr<ThermoParams> thermo_params_;
  std::shared_ptr<Criteria> criteria_;
  TrialFactory extra_trials_;
  std::shared_ptr<ProgressReport> report_;
};

inline std::shared_ptr<SeekNumParticles> MakeSeekNumParticles(
    const int num,
    const argtype &args = argtype()) {
  return std::make_shared<SeekNumParticles>(num, args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_
