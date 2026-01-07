
#ifndef FEASST_MORPH_TRIAL_MORPH_H_
#define FEASST_MORPH_TRIAL_MORPH_H_

#include <vector>
#include <string>
#include <memory>
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Morph random particle(s) of given type into a different particle type(s).
  Typically requires the use of a reference index if multiple particles are to
  be morphed simultaneously.
  See ComputeMorph.
 */
class TrialMorph : public TrialFactoryNamed {
 public:
  //@{
  /** @name Arguments
    - particle_type: name of particle type that will be morphed.
      Multiple types can be given with comma-separated.
    - particle_type_morph: name of particle type to morph into.
      Multiple arguments must correspond with the above.
      For example, "particle_type 1,2 particle_type_morph 3,4"
      changes 1->2 and 3->4.
    - Trial arguments.
    - TrialSelectParticle arguments (but particle_type is specified as above).
   */
  explicit TrialMorph(argtype args = argtype());
  explicit TrialMorph(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialMorph>(args); }
  virtual ~TrialMorph() {}
  //@}
};

class TrialMorphOneWay : public Trial {
 public:
  explicit TrialMorphOneWay(argtype args = argtype());
  explicit TrialMorphOneWay(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialMorphOneWay>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialMorphOneWay>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialMorphOneWay(std::istream& istr);
  virtual ~TrialMorphOneWay() {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_MORPH_TRIAL_MORPH_H_
