
#ifndef FEASST_MONTE_CARLO_COMPUTE_AVB2_H_
#define FEASST_MONTE_CARLO_COMPUTE_AVB2_H_

#include <memory>
#include <vector>
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 */
class ComputeAVB2 : public TrialComputeMove {
 public:
  /**
    args:
    - probability_out_to_in: probability to attempt an out->in move
      (default: 0.5).
    - out_to_in: true if out->in move, otherwise false.
   */
  explicit ComputeAVB2(argtype args = argtype());
  explicit ComputeAVB2(argtype * args);

  /// Return the probability of an out to in move.
  double probability_out_to_in() const { return p_bias_; }

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeAVB2(std::istream& istr);
  virtual ~ComputeAVB2() {}

 protected:
  void serialize_compute_avb2_(std::ostream& ostr) const;
  double p_bias_;
  bool out_to_in_;
};

inline std::shared_ptr<ComputeAVB2> MakeComputeAVB2(argtype args = argtype()) {
  return std::make_shared<ComputeAVB2>(args);
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_COMPUTE_AVB_H_
