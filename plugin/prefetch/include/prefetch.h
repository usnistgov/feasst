
#ifndef FEASST_PREFETCH_PREFETCH_H_
#define FEASST_PREFETCH_PREFETCH_H_

#include <string>
#include <vector>
#include <memory>
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Thread;

// HWH: consider a parallel while:
// https://cvw.cac.cornell.edu/OpenMP/whileloop

/**
  Define a pool of threads, each with their own MonteCarlo object and
  relevant quantities for reversion.
*/
class Pool {
 public:
  void set_index(const int index) {
    index_ = index; }
  int index() const { return index_; }
  void set_ln_prob(const double ln_prob) { ln_prob_ = ln_prob; }
  double ln_prob() const { return ln_prob_; }
  void set_accepted(const bool accepted) { accepted_ = accepted; }
  bool accepted() const { return accepted_; }
  void set_auto_rejected(const bool auto_rejected) { auto_rejected_ = auto_rejected; }
  bool auto_rejected() const { return auto_rejected_; }
  void set_endpoint(const bool endpoint) { endpoint_ = endpoint; }
  bool endpoint() const { return endpoint_; }
  const std::string str() const;
  std::unique_ptr<MonteCarlo> mc;

 private:
  int index_;
  double ln_prob_;
  bool accepted_;
  bool auto_rejected_ = false;
  bool endpoint_ = true;
};

/**
  Farm a trial to each processor, then reconstruct the serial Markov chain.

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
  /**
    args:
    - trials_per_check: number of steps between check (default: 1e6)
    - load_balance: batches contain all of the same trial type (default: false).
      This violates detailed balance, and is known in cases of high acceptance
      to give erroneous results.
      Only use load_balance for equilibration and never for production simulations.
    - synchronize: synchronize data with accepted thread (default: false).
    - ghost: update transition matrix even for trials after acceptance (default: false).
   */
  explicit Prefetch(argtype args = argtype());

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
  void run_num_trials(int num_trials) override;
  void run_until_num_particles(const int num_particles,
                               const std::string& particle_type,
                               const int configuration_index) override;
  void run_until_file_exists(const std::string& file_name,
                             const int trials_per_file_check) override;
  void run_until_complete() override;

  void serialize(std::ostream& ostr) const override;
  explicit Prefetch(std::istream& istr);
  virtual ~Prefetch() {}

 protected:
  void attempt_(int num_trials, TrialFactory * trial_factory, Random * random) override;
  void run_until_complete_(TrialFactory * trial_factory, Random * random) override;

 private:
  bool is_activated_;
  bool is_synchronize_;
  int trials_per_check_;
  int trials_since_check_ = 0;
  bool load_balance_;
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
