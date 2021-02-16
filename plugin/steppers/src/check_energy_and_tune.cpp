#include "utils/include/serialize.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy_and_tune.h"

namespace feasst {

CheckEnergyAndTune::CheckEnergyAndTune(argtype args) : ModifyFactory() {
  argtype check_args = args;
  add(MakeCheckEnergy(check_args));

  // remove tolerance and apply the same args to Tune
  str("tolerance", &args, "?");
  add(MakeTuner(args));
}

}  // namespace feasst
