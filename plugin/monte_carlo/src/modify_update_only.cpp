#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/modify_update_only.h"

namespace feasst {

ModifyUpdateOnly::ModifyUpdateOnly(argtype * args) : Modify(args) {
  // disable write
  Modify::set_trials_per_write(-1);

  // parse
  if (used("trials_per", *args)) {
    WARN("ModifyUpdateOnly::trials_per is deprecated. Use trials_per_update.");
    set_trials_per(integer("trials_per", args));
  }
  ASSERT(output_file().empty(),
    "ModifyUpdateOnly does not use the argument output_file.");
}

void ModifyUpdateOnly::set_trials_per_write(const int trials) {
  ERROR("This modify is update only.");
}

}  // namespace feasst
