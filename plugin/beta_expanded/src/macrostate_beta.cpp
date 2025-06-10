#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/histogram.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "beta_expanded/include/macrostate_beta.h"

namespace feasst {

MacrostateBeta::MacrostateBeta(const Histogram& histogram,
    argtype * args) : Macrostate(histogram, args) {
  class_name_ = "MacrostateBeta";
}
MacrostateBeta::MacrostateBeta(const Histogram& histogram,
    argtype args) : Macrostate(histogram, args) {
  feasst_check_all_used(args);
}
MacrostateBeta::MacrostateBeta(argtype args) :
    MacrostateBeta(Histogram(&args), &args) {
  feasst_check_all_used(args);
}
MacrostateBeta::MacrostateBeta(argtype * args) :
  MacrostateBeta(Histogram(args), args) {}

FEASST_MAPPER(MacrostateBeta,);

std::shared_ptr<Macrostate> MacrostateBeta::create(std::istream& istr) const {
  return std::make_shared<MacrostateBeta>(istr);
}

std::shared_ptr<Macrostate> MacrostateBeta::create(argtype * args) const {
  return std::make_shared<MacrostateBeta>(args);
}

MacrostateBeta::MacrostateBeta(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1048, "version mismatch: " << version);
}

void MacrostateBeta::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(1048, ostr);
}

double MacrostateBeta::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) {
  return system.thermo_params().beta();
}

}  // namespace feasst
