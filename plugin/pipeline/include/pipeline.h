
#ifndef FEASST_PIPELINE_PIPELINE_H_
#define FEASST_PIPELINE_PIPELINE_H_

#include <string>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

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
class Pipeline : public MonteCarlo {
 public:
  Pipeline() {
    activate_pipeline();
  }

  void activate_pipeline(const bool active = true) { is_activated_ = active; }

  void serialize(std::ostream& ostr) const override { ERROR("not implemented"); }
  Pipeline(std::istream& istr) : MonteCarlo(istr) {}

  virtual ~Pipeline() {}

 protected:
  void attempt_(const int num_trials, TrialFactory * trial_factory, Random * random) override;

 private:
  // void distribute_for_();

  // Distribute pool of trial(s) to each processor
  // use a omp parallel while constructed as described in
  // https://cvw.cac.cornell.edu/OpenMP/whileloop
//  void distribute_while_();

  int num_threads_;
  bool is_activated_ = true;
  int steps_per_check_ = 1e4;

  std::vector<Pool> pool_;

  // Pick a clone based on jthread. If i==j, use this object instead.
  MonteCarlo * clone_(const int ithread, const int jthread);
};

}  // namespace feasst

#endif  // FEASST_PIPELINE_PIPELINE_H_
