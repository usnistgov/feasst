
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

// HWH
/**
  Defines the macrostate to be the total potential energy of the system.
  WARNING: this class has not been tested at all.
 */
class MacrostateEnergy : public Macrostate {
 public:
  //@{
  /** @name Arguments
    - particle_type: number of particles of type. If -1 (default), count all
      types.
    - Macrostate arguments.
  */
  explicit MacrostateEnergy(argtype args = argtype());
  explicit MacrostateEnergy(argtype * args) :
    MacrostateEnergy(Histogram(args), args) {}

  //@}
  /** @name Public Functions
   */
  //@{

  /// Arguments as described above, but with explicit histogram object.
  MacrostateEnergy(const Histogram& histogram, argtype args = argtype());
  MacrostateEnergy(const Histogram& histogram, argtype * args);

  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateEnergy(std::istream& istr);
  virtual ~MacrostateEnergy() {}
  //@}
};

inline std::shared_ptr<MacrostateEnergy> MakeMacrostateEnergy(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostateEnergy>(histogram, args);
}

inline std::shared_ptr<MacrostateEnergy> MakeMacrostateEnergy(
    argtype args = argtype()) {
  return std::make_shared<MacrostateEnergy>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_
