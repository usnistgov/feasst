
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_

#include "monte_carlo/include/constrain_num_particles.h"
#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total number of particles in the system.
 */
class MacrostateNumParticles : public Macrostate {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - particle_type: number of particles of type. If -1 (default), count all
      types.
    - Macrostate arguments.
  */
  explicit MacrostateNumParticles(argtype args = argtype());
  explicit MacrostateNumParticles(argtype * args) :
    MacrostateNumParticles(Histogram(args), args) {}

  //@}
  /** @name Public Functions
   */
  //@{

  /// Arguments as described above, but with explicit histogram object.
  MacrostateNumParticles(const Histogram& histogram, argtype args = argtype());
  MacrostateNumParticles(const Histogram& histogram, argtype * args);

  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateNumParticles(std::istream& istr);
  virtual ~MacrostateNumParticles() {}

  //@}
 private:
  // int particle_type_;
  ConstrainNumParticles num_;
};

inline std::shared_ptr<MacrostateNumParticles> MakeMacrostateNumParticles(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostateNumParticles>(histogram, args);
}

inline std::shared_ptr<MacrostateNumParticles> MakeMacrostateNumParticles(
    argtype args = argtype()) {
  return std::make_shared<MacrostateNumParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
