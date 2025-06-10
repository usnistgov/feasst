
#ifndef FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
#define FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_

#include <memory>
#include <string>
#include <map>
#include "flat_histogram/include/macrostate.h"

namespace feasst {

class Acceptance;
class Criteria;
class System;

typedef std::map<std::string, std::string> argtype;

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
  explicit MacrostateBeta(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{
  MacrostateBeta(const Histogram& histogram, argtype args = argtype());
  MacrostateBeta(const Histogram& histogram, argtype * args);
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateBeta(std::istream& istr);
  virtual ~MacrostateBeta() {}
  //@}
};

}  // namespace feasst

#endif  // FEASST_BETA_EXPANDED_MACROSTATE_BETA_H_
