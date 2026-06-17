#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "monte_carlo/include/analyze_update_only.h"

namespace feasst {

AnalyzeUpdateOnly::AnalyzeUpdateOnly(argtype * args) : Analyze(args) {
  ASSERT(trials_per_write() == 1,
    "AnalyzeUpdateOnly does not use the argument trials_per_write.");
  ASSERT(output_file().empty(),
    "AnalyzeUpdateOnly does not use the argument output_file.");

  // disable write
  Analyze::set_trials_per_write(-1);

  // parse
  if (used("trials_per", *args)) {
    set_trials_per(integer("trials_per", args));
  }
}

void AnalyzeUpdateOnly::set_trials_per_write(const int trials) {
  ERROR("This analyze is update only.");
}

}  // namespace feasst
