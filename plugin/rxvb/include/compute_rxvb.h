
#ifndef FEASST_RXVB_COMPUTE_RXVB_H_
#define FEASST_RXVB_COMPUTE_RXVB_H_

#include <memory>
#include <vector>
#include "configuration/include/select.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 */
class ComputeRXVB : public TrialComputeMove {
 public:
  /**
    args:
    - probability_avb: probability to attempt an AVB rxn (default: 0.5).
    - avb: true if AVB rxn, otherwise reverse rxn without AVB (default: true).
    - target_particle_type
    - particle_type
    - target_particle_type_morph
    - particle_type_morph
   */
  explicit ComputeRXVB(argtype args = argtype());
  explicit ComputeRXVB(argtype * args);

  /// Return the probability of an AVB rxn.
  double probability_avb() const { return p_bias_; }

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override;

  // serialize
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ComputeRXVB(std::istream& istr);
  virtual ~ComputeRXVB() {}

 protected:
  void serialize_compute_rxnavb_(std::ostream& ostr) const;
  double p_bias_;
  bool avb_;
  int target_particle_type_ = -2;
  int particle_type_ = -2;
  int target_particle_type_morph_ = -2;
  int particle_type_morph_ = -2;
  std::string target_particle_type_name_;
  std::string particle_type_name_;
  std::string target_particle_type_morph_name_;
  std::string particle_type_morph_name_;

  // temporary and not serialized
  Select neighbors_;
  int num_particles_of_type_(const TrialSelect& select, const Configuration& config);
};

}  // namespace feasst

#endif  // FEASST_RXVB_COMPUTE_RXVB_H_
