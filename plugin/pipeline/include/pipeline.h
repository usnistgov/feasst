
#ifndef FEASST_PIPELINE_PIPELINE_H_
#define FEASST_PIPELINE_PIPELINE_H_

#include <string>
#include <vector>
#include <memory>
#include "utils/include/arguments.h"
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
  Pipeline(
    /**
      steps_per_check : number of steps between check (default: 1e5)
     */
    const argtype& args = argtype()) {
    activate_pipeline();
    Arguments args_(args);
    steps_per_check_ = args_.key("steps_per_check").dflt("100000").integer();
  }

  void activate_pipeline(const bool active = true) { is_activated_ = active; }
  void reset_trial_stats() override {
    MonteCarlo::reset_trial_stats();
    for (Pool& pool : pool_) {
      pool.mc.reset_trial_stats();
    }
  }

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
  bool is_activated_;
  int steps_per_check_;

  std::vector<Pool> pool_;

  // Pick a clone based on jthread. If i==j, use this object instead.
  MonteCarlo * clone_(const int ithread, const int jthread);
};

inline std::shared_ptr<Pipeline> MakePipeline(const argtype& args = argtype()) {
  return std::make_shared<Pipeline>(args);
}

}  // namespace feasst

#endif  // FEASST_PIPELINE_PIPELINE_H_
