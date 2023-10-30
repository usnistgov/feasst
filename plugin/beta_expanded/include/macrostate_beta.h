
#ifndef FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
#define FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the inverse temperature, \f$\beta\f$.
 */
class MacrostateBeta : public Macrostate {
 public:
  //@{
  /** @name Arguments
    - Macrostate arguments.
   */
  explicit MacrostateBeta(argtype args = argtype());
  explicit MacrostateBeta(argtype * args) :
    MacrostateBeta(Histogram(args), args) {}
  //@}
  /** @name Public Functions
   */
  //@{
  MacrostateBeta(const Histogram& histogram, argtype args = argtype());
  MacrostateBeta(const Histogram& histogram, argtype * args);
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override {
    return system.thermo_params().beta(); }
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateBeta(std::istream& istr);
  virtual ~MacrostateBeta() {}
  //@}
};

inline std::shared_ptr<MacrostateBeta> MakeMacrostateBeta(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostateBeta>(histogram, args);
}

inline std::shared_ptr<MacrostateBeta> MakeMacrostateBeta(
    argtype args = argtype()) {
  return std::make_shared<MacrostateBeta>(args);
}

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
