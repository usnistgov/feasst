
#ifndef FEASST_GROWTH_EXPANDED_MACROSTATE_GROWTH_EXPANDED_H_
#define FEASST_GROWTH_EXPANDED_MACROSTATE_GROWTH_EXPANDED_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total number of particles in the system
  plus the partial fraction of the growing particle.
 */
class MacrostateGrowthExpanded : public Macrostate {
 public:
  MacrostateGrowthExpanded(const Histogram& histogram,
    const argtype& args = argtype()) : Macrostate(histogram, args) {}
  double value(const System* system,
               const Criteria* criteria,
               const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateGrowthExpanded(std::istream& istr);
  virtual ~MacrostateGrowthExpanded() {}

 private:
  const std::string class_name_ = "MacrostateGrowthExpanded";
};

inline std::shared_ptr<MacrostateGrowthExpanded> MakeMacrostateGrowthExpanded(
    const Histogram& histogram, const argtype& args = argtype()) {
  return std::make_shared<MacrostateGrowthExpanded>(histogram, args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_MACROSTATE_GROWTH_EXPANDED_H_
