
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total potential energy of the system.
 */
class MacrostateEnergy : public Macrostate {
 public:
  /**
   args:
   - particle_type: number of particles of type. If -1 (default), count all
     types.
  */
  MacrostateEnergy(const Histogram& histogram, const argtype& args = argtype());
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateEnergy(std::istream& istr);
  virtual ~MacrostateEnergy() {}

 private:
  const std::string class_name_ = "MacrostateEnergy";
};

inline std::shared_ptr<MacrostateEnergy> MakeMacrostateEnergy(
    const Histogram& histogram, const argtype& args = argtype()) {
  return std::make_shared<MacrostateEnergy>(histogram, args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_ENERGY_H_
