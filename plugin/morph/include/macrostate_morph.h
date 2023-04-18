
#ifndef FEASST_GROWTH_EXPANDED_MACROSTATE_MORPH_H_
#define FEASST_GROWTH_EXPANDED_MACROSTATE_MORPH_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "flat_histogram/include/macrostate.h"

namespace feasst {

// HWH if type is -1, spaceholder, don't do anything
// HWH but everything after -1 must be more -1 within the same step
// HWH implementation was stalled because macrostate doesn't distinguish
// HWH between different steps
/**
  For use with TrialMorphExpanded.

  The particle_type_growth_sequence (PTGS) argument is as follows.

  The first index represents individual macrostates.
  Thus, if the first index is of size 2, the Histogram width should be set equal
  to 1/PTGS, which means that macrostates would be 0, 0.5, 1.0, 1.5, ...

  The second index represents the (multiple) particles that are to be morphed.

  For example, PTGS={{1},{0}}, would mean to insert a particle of type 1 to
  transition to a half macrostate value, and morph particles of type 1 into
  type 0 to obtain integer macrostates.
  The reverse at integer macrostates would be to randomly select a particle
  of type 0 and shrink it into a particle of type 1.
  The reverse at half macrostate would be to delete particles of type 1.

  Thus, assuming particle type 1 is non-physical and introduced only to help
  sampling, the integer macrostates should contain only particles of type 0.
  Thus, there should only be particles of type 1 present at half macrostates.

  In a more complex example for divalent ions, PTGS={{2, 3, 3}, {0, 1, 1}} if
  particles of type 0 have twice the charge of particles of type 1.
  In this case, we begin by adding partial particles of types 2 and 3 at the
  half macrostate.
  Those particles are then morphed into full particles of types 1 and 0 when
  transitioning to integer macrostate values.
 */
class MacrostateMorph : public Macrostate {
 public:
  MacrostateMorph(
    const std::vector<std::vector<int> > particle_type_growth_sequence,
    const Histogram& histogram,
    const argtype& args = argtype());
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit MacrostateMorph(std::istream& istr);
  virtual ~MacrostateMorph() {}

 private:
  std::vector<std::vector<int> > grow_seq_;
  std::vector<ConstrainNumParticles> num_first_;
  bool is_row_all_populated_(const int row, const Matrix& matrix) const;
};

inline std::shared_ptr<MacrostateMorph> MakeMacrostateMorph(
    const std::vector<std::vector<int> > particle_type_growth_sequence,
    const Histogram& histogram, const argtype& args = argtype()) {
  return std::make_shared<MacrostateMorph>(particle_type_growth_sequence,
    histogram, args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_MACROSTATE_MORPH_H_
