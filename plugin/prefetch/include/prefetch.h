
#ifndef FEASST_PREFETCH_PREFETCH_H_
#define FEASST_PREFETCH_PREFETCH_H_

#include <string>
#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

// HWH: consider a parallel while:
// https://cvw.cac.cornell.edu/OpenMP/whileloop

/**
  Define a pool of threads, each with their own MonteCarlo object and
  relevant quantities for reversion.
*/
class Pool {
 public:
  void set_index(const int index) {
    DEBUG("index " << index);
    index_ = index; }
  int index() const { return index_; }
  void set_ln_prob(const double ln_prob) { ln_prob_ = ln_prob; }
  double ln_prob() const { return ln_prob_; }
  void set_accepted(const bool accepted) { accepted_ = accepted; }
  bool accepted() const { return accepted_; }

  const std::string str() const {
    std::stringstream ss;
    ss << index_ << " " << ln_prob_ << " " << accepted_;
    return ss.str();
  }

  MonteCarlo mc;

 private:
  int index_;
  double ln_prob_;
  bool accepted_;
};

/**
  Farm a trial to each processor, then reconstruct the serial Markov chain.
 */
class Prefetch : public MonteCarlo {
 public:
  /**
    args:
    - steps_per_check: number of steps between check (default: 1e5)
    - load_balance: batches contain all of the same trial type (default: true).
      This violates detailed balance, but "local detailed balance" is still
      satisfied.
      https://doi.org/10.1063/1.477973
   */
  explicit Prefetch(const argtype& args = argtype());

  /// Return the number of steps between checking equality of threads.
  int steps_per_check() const { return steps_per_check_; }

  /// Activate prefetch.
  void activate_prefetch(const bool active = true) { is_activated_ = active; }

  /// Reset stats of trials of all threads.
  void reset_trial_stats() override;

  // public interface for unit testing only
  const std::vector<Pool>& pool() const { return pool_; }
  // Pick a clone based on ithread.
  // If ithread == 0, return self. Otherwise, return pool_.
  MonteCarlo * clone_(const int ithread);

  void serialize(std::ostream& ostr) const override;
  explicit Prefetch(std::istream& istr);
  virtual ~Prefetch() {}

 protected:
  void attempt_(int num_trials, TrialFactory * trial_factory, Random * random) override;
  void run_until_complete_(TrialFactory * trial_factory, Random * random) override;

 private:
  bool is_activated_;
  int steps_per_check_;
  int steps_since_check_ = 0;
  bool load_balance_;

  // temporary
  int num_threads_;
  std::vector<Pool> pool_;

  void create(std::vector<Pool> * pool);
};

inline std::shared_ptr<Prefetch> MakePrefetch(const argtype& args = argtype()) {
  return std::make_shared<Prefetch>(args);
}

}  // namespace feasst

#endif  // FEASST_PREFETCH_PREFETCH_H_
