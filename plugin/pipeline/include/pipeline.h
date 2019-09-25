
#ifndef FEASST_PIPELINE_PIPELINE_H_
#define FEASST_PIPELINE_PIPELINE_H_

#include <string>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

class Pool {
 public:
  void set_index(const int index) {
    INFO("index " << index);
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

  private:
   int index_;
   double ln_prob_;
   bool accepted_;
};

/**
  A trial attempt farms trials to processors then reconstructs the serial
  Markov chain.
 */
class Pipeline : public MonteCarlo {
 public:
  Pipeline() {}
//    set_num_pool();
//  }

  /// Set number of trials in pool. If -1, set to number of threads.
//  void set_num_pool(const int num = -1) { num_pool_ = num; }

  Pipeline(std::istream& istr) : MonteCarlo(istr) {}

 private:
  // void distribute_for_();

  // Distribute pool of trial(s) to each processor
  // use a omp parallel while constructed as described in
  // https://cvw.cac.cornell.edu/OpenMP/whileloop
//  void distribute_while_();

  std::vector<MonteCarlo> clones_;
//  int num_pool_;
  int num_threads_;

  std::vector<Pool> trial_pool_;

  void attempt_(
      System * system,
      Criteria * criteria,
      TrialFactory * trial_factory,
      AnalyzeFactory * analyze_factory,
      ModifyFactory * modify_factory,
      Checkpoint * checkpoint,
      Random * random) override;

  // Pick a clone based on jthread. If i==j, use this object instead.
  MonteCarlo * clone_(const int ithread, const int jthread);
};

}  // namespace feasst

#endif  // FEASST_PIPELINE_PIPELINE_H_
