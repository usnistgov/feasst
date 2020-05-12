#include "steppers/include/utils.h"
#include "utils/include/utils_io.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

void add_common_steppers(MonteCarlo * monte_carlo, const argtype& args) {
  Arguments args_(args);
  const std::string steps_per = args_.key("steps_per").dflt(str(1e6)).str();
  const std::string file_app = args_.key("file_append").str();
  monte_carlo->add(MakeLog(
   {{"steps_per", steps_per},
    {"file_name", file_app + "_log.txt"}}));
  monte_carlo->add(MakeMovie(
   {{"steps_per", steps_per},
    {"file_name", file_app + "_movie.xyz"}}));
  monte_carlo->add(MakeCheckEnergy(
   {{"steps_per", steps_per},
    {"tolerance", args_.key("tolerance").dflt(str(1e-10)).str()}}));
  monte_carlo->add(MakeTuner({{"steps_per", str(steps_per)}}));
}

}  // namespace feasst
