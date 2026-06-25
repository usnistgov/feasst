
#ifndef FEASST_PREFETCH_PREFETCH_H_
#define FEASST_PREFETCH_PREFETCH_H_

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include "monte_carlo/include/monte_carlo.h"
#include "prefetch/include/pool.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Set the number of threads using the BASH environmental command:

  export OMP_NUM_THREADS=2

  Prefetch farma a trial to each processor, then reconstruct the serial Markov chain.

  Although its quite a simple strategy, the trickiest part is efficiently
  reproducing the accepted perturbation on all threads (e.g., the synchronizing
  step) in a way that's general to all MC trials (and potential functions/
  optimizations).
  This can get really complicated by cell lists, neighbor lists and k-space.
  For the synchronize step, I've used two different strategies, sometimes both
  at the same time, which span the range of ease of implementation vs efficiency.

  1. Cache every random number and energy term generated during the trial, and
  if that trial is the accepted one, revert the others, then feed those cached
  RNG/energies back to the other threads to reproduce the same perturbations.
  There are some issues with FEASST here because the energy calcs are also part
  of the neighbor/cell/kvector updates.
  Also, this means going through each config bias step, etc, so its not the most
  efficient synchronization possible.

  2. A more efficient synchronization method is to give objects a synchronize
  method which takes as input a reference to the base class its trying to copy
  and the list of particles/sites that were changed in the trial.
  For example, in Ewald, VisitModel is the base class of Ewald, and Select has
  the list of particles/sites that were perturbed.
  The base class also needs a data structure (SynchronizeData) that all the
  derived classes use.
  In the most complex cases, the data is separated between ones that are
  automatically copied every time, or ones manually choosen based on which
  sites were perturbed.

  Prefetch is not used for the until_num_particles argument in Run.
 */
class Prefetch : public MonteCarlo {
 public:
  //@{
  /** @name Arguments
    - trials_per_check: number of steps between check (default: 1e6)
    - load_balance: number of trials in a batch of attempted trials of the same
      trial type (default: -1).
      If < 1, do not load balance.
      If == 1 (not recommended), load balance more efficiently by breaking
      detailed balance, where each batch is simply the number of processors
      regardless of first accepted trial.
      If > 1, load balance by performing this many attempted trials, not
      counting trials after the first accepted.
    - synchronize: synchronize data with accepted thread (default: false).
    - ghost: update transition matrix even for trials after acceptance (default: false).
   */
  explicit Prefetch(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{


  /// Return the number of steps between checking equality of threads.
  int trials_per_check() const { return trials_per_check_; }

  /// Activate prefetch.
  void activate_prefetch(const bool active = true) { is_activated_ = active; }

  /// Reset stats of trials of all threads.
  void reset_trial_stats() override;

  // public interface for unit testing only
  const std::vector<Pool>& pool() const { return pool_; }
  // Pick a clone based on ithread.
  // If ithread == 0, return self. Otherwise, return pool_.
  MonteCarlo * clone_(const int ithread);

  /// Perform an Action on all processors.
  void run(std::shared_ptr<Action> action) override;

  /// Run a number of trials.
  void run_num_trials(int64_t num_trials) override;
  void run_until_num_particles(const int num_particles,
                               const std::string& particle_type,
                               const int configuration_index) override;
  void run_until_file_exists(const std::string& file_name,
                             const int trials_per_file_check) override;
  void run_until_volume(const double volume, const int configuration_index) override;
  void run_until_complete() override;

  void serialize(std::ostream& ostr) const override;
  explicit Prefetch(std::istream& istr);
  virtual ~Prefetch() {}

  //@}
 protected:
  void attempt_(int num_trials, TrialFactory * trial_factory, Random * random) override;
  void run_until_complete_(TrialFactory * trial_factory, Random * random) override;

 private:
  bool is_activated_;
  bool is_synchronize_;
  int trials_per_check_;
  int trials_since_check_ = 0;
  int trials_since_balance_ = 0;
  int load_balance_;
  int load_balance_trial_;
  bool ghost_;

  // temporary
  int num_threads_;
  std::vector<Pool> pool_;

  void create(std::vector<Pool> * pool);
};

inline std::shared_ptr<Prefetch> MakePrefetch(argtype args = argtype()) {
  return std::make_shared<Prefetch>(args);
}

}  // namespace feasst

#endif  // FEASST_PREFETCH_PREFETCH_H_
