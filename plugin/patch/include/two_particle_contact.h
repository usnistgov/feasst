
#ifndef FEASST_PATCH_TWO_PARTICLE_CONTACT_H_
#define FEASST_PATCH_TWO_PARTICLE_CONTACT_H_

#include <vector>
#include <memory>
#include "math/include/formula.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/action.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class TwoParticleContact;

class TwoParticleContactObjective : public Formula {
 public:
  TwoParticleContactObjective(TwoParticleContact * two_particle_contact, System * system);
  double evaluate(const double distance) const override;
 private:
  TwoParticleContact * two_particle_contact_;
  System * system_;
};

/**
  Generate the contact distance between two particles.
  The current implementation assumes these are hard particles.
  Thus, the energy can only be zero or NEAR_INFINITY.
 */
class TwoParticleContact : public Action {
 public:
  //@{
  /** @name Arguments
    - dimension: dimension to move particles along to find contact (default: 0).
    - output_file: name of file to output contact distance.
      Do not print if empty (default: empty).
    - tolerance: tolerance of minimum distance search (default: 1e-6).
   */
  explicit TwoParticleContact(argtype args = argtype());
  explicit TwoParticleContact(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy(const double distance, System * system);
  void update_xyz(const double distance, System * system);
  void revert(System * system);
  double contact_distance() const { return contact_distance_; }
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<TwoParticleContact>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<TwoParticleContact>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TwoParticleContact(std::istream& istr);
  virtual ~TwoParticleContact() {}

  //@}
 private:
  std::string output_file_;
  int dimension_;
  double tolerance_;

  // no need to serialize
  double contact_distance_;
  std::shared_ptr<PerturbTranslate> translate_;
  std::shared_ptr<TrialSelectParticle> select_;
  Position com1_;
};

}  // namespace feasst

#endif  // FEASST_PATCH_TWO_PARTICLE_CONTACT_H_
