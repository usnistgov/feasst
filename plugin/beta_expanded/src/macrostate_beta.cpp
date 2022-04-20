#include "utils/include/serialize.h"
#include "beta_expanded/include/macrostate_beta.h"

namespace feasst {

MacrostateBeta::MacrostateBeta(const Histogram& histogram,
    argtype * args) : Macrostate(histogram, args) {
  class_name_ = "MacrostateBeta";
}
MacrostateBeta::MacrostateBeta(const Histogram& histogram,
    argtype args) : Macrostate(histogram, args) {
  check_all_used(args);
}
MacrostateBeta::MacrostateBeta(argtype args) :
    MacrostateBeta(Histogram(&args), &args) {
  check_all_used(args);
}

class MapMacrostateBeta {
 public:
  MapMacrostateBeta() {
    auto hist = MakeHistogram({{"width", "1"}, {"max", "1"}});
    MacrostateBeta(*hist).deserialize_map()["MacrostateBeta"] =
      MakeMacrostateBeta(*hist);
  }
};

static MapMacrostateBeta mapper_ = MapMacrostateBeta();

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

}  // namespace feasst
