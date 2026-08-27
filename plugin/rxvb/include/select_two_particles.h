
#ifndef FEASST_RXVB_SELECT_TWO_PARTICLES_H_
#define FEASST_RXVB_SELECT_TWO_PARTICLES_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 * Randomly Select two particles of different types.
 */
class SelectTwoParticles : public TrialSelect {
 public:
  //@{
  /** @name Arguments
    - particle_types: two comma-separated names of particles types that must be
      different.
   */
  explicit SelectTwoParticles(argtype args = argtype());
  explicit SelectTwoParticles(argtype * arg);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(System * system) override;
  bool select(const Select& perturbed,
    System* system,
    Random * random,
    TrialSelect * previous_select) override;
  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectTwoParticles(std::istream& istr);
  virtual ~SelectTwoParticles();

  //@}
 private:
  std::vector<int> pts_;
  std::vector<std::string> pt_names_;
  std::vector<int> group_indices_;
};

}  // namespace feasst

#endif  // FEASST_RXVB_SELECT_TWO_PARTICLES_H_
