
#ifndef FEASST_MONTE_CARLO_CONSTRAIN_NUM_PARTICLES_H_
#define FEASST_MONTE_CARLO_CONSTRAIN_NUM_PARTICLES_H_

#include <memory>
#include "monte_carlo/include/constraint.h"

namespace feasst {

/**
  Constrain the number of the number of particles.
 */
class ConstrainNumParticles : public Constraint {
 public:
  //@{
  /** @name Arguments
    - maximum: maximum number of particles. If -1, no limit (default: -1).
    - minimum: minimum number of particles (default: 0).
    - type: particle type name. If empty (default), all types.
    - configuration_index: index of configuration (default: 0).
   */
  explicit ConstrainNumParticles(argtype args = argtype());
  explicit ConstrainNumParticles(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the particle type.
  int type() const { return type_; }

  bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) override;

  /// Return the number of particles, taking into account the potential to shift
  /// (e.g., add or delete) from a trial move.
  int num_particles(const System& system, const Acceptance& acceptance);

  std::shared_ptr<Constraint> create(std::istream& istr) const override {
    return std::make_shared<ConstrainNumParticles>(istr); }
  std::shared_ptr<Constraint> create(argtype * args) const override {
    return std::make_shared<ConstrainNumParticles>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ConstrainNumParticles(std::istream& istr);
  virtual ~ConstrainNumParticles() {}

  //@}
 protected:
  void serialize_constrain_num_particles_(std::ostream& ostr) const;
  int maximum_;
  int minimum_;
  int type_;
  std::string type_name_;
  int configuration_index_;
};

inline std::shared_ptr<ConstrainNumParticles> MakeConstrainNumParticles(
    argtype args = argtype()) {
  return std::make_shared<ConstrainNumParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CONSTRAIN_NUM_PARTICLES_H_
