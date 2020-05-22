
#ifndef FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_
#define FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_

#include <memory>
#include "utils/include/arguments.h"

namespace feasst {

class Criteria;
class TrialAdd;
class MonteCarlo;

/**
  Initialize the number of particles in MonteCarlo.
  Use the existing trials, but reset their stats once the goal is reached.
 */
class SeekNumParticles {
 public:
  /**
    args:
    - particle_type: type of particle to add (default: 0)
    - max_attempts: maximum number of trial attempts for seek (default: 1e8).
   */
  SeekNumParticles(const int num, const argtype& args = argtype());

  /// Optionally use a temporary Metropolis criteria with given args.
  SeekNumParticles with_metropolis(const argtype& args);

  /// Use an extra and temporary TrialAdd between every existing Trial.
  /// If particle_type is not included in args, assume 0.
  SeekNumParticles with_trial_add(const argtype& args = argtype());

  /// Perform the seek.
  void run(MonteCarlo * monte_carlo);

 private:
  const int num_;
  int max_attempts_;
  int particle_type_;
  std::shared_ptr<Criteria> criteria_;
  std::shared_ptr<TrialAdd> add_;
};

inline std::shared_ptr<SeekNumParticles> MakeSeekNumParticles(
    const int num,
    const argtype &args = argtype()) {
  return std::make_shared<SeekNumParticles>(num, args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_SEEK_NUM_PARTICLES_H_
