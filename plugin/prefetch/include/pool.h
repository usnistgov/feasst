#ifndef FEASST_PREFETCH_POOL_H_
#define FEASST_PREFETCH_POOL_H_

#include <memory>
#include <string>
#include <sstream>

namespace feasst {

class MonteCarlo;

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

}  // namespace feasst

#endif  // FEASST_PREFETCH_POOL_H_
