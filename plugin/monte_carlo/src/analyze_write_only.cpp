#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "monte_carlo/include/analyze_write_only.h"

namespace feasst {

AnalyzeWriteOnly::AnalyzeWriteOnly(argtype * args) : Analyze(args) {
  // disable update
  Stepper::set_trials_per_update(-1);

  // parse
  if (used("trials_per", *args)) {
    WARN("AnalyzeWriteOnly::trials_per is deprecated. Use trials_per_write.");
    set_trials_per(integer("trials_per", args));
  }
}

void AnalyzeWriteOnly::set_trials_per_update(const int trials) {
  ERROR("This analyze is write only.");
}

}  // namespace feasst
