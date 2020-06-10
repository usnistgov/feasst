#include "utils/include/serialize.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy_and_tune.h"

namespace feasst {

CheckEnergyAndTune::CheckEnergyAndTune(const argtype& args) : ModifyFactory() {
  add(MakeCheckEnergy(args));
  add(MakeTuner(Arguments().remove("tolerance", args)));
}

}  // namespace feasst
