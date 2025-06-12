#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/file_vmd.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/movie.h"

namespace feasst {

FEASST_MAPPER(Movie, argtype({{"output_file", "place_holder"}}));

Movie::Movie(argtype * args) : AnalyzeWriteOnly(args) {
  set_append();
  ASSERT(!output_file().empty(), "output_file argument is required");
  args->insert({"append", "true"}); // always append
  xyz_ = std::make_unique<FileXYZ>(args);
  vmd_ = std::make_unique<FileVMD>(args);
}
Movie::Movie(argtype args) : Movie(&args) { feasst_check_all_used(args); }
Movie::~Movie() {}

void Movie::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  const System& system = mc->system();
  const std::string name = output_file(mc->criteria());
  ASSERT(!name.empty(), "output_file argument is required. Did you forget to " <<
    "Analyze::set_output_file()?");

  // write xyz
  if (state() == mc->criteria().state()) {
    xyz_->write(name, configuration(system));
  }

  // write vmd
  std::stringstream ss;
  ss << name << ".vmd";
  vmd_->write(ss.str(), configuration(system), name);
}

std::string Movie::write(const MonteCarlo& mc) {
  // ensure the following order matches the header from initialization.
  xyz_->write(output_file(mc.criteria()), configuration(mc.system()));
  return std::string("");
}

void Movie::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(536, ostr);
  feasst_serialize(xyz_, ostr);
  feasst_serialize(vmd_, ostr);
}

Movie::Movie(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 536, "version mismatch:" << version);
  feasst_deserialize(xyz_, istr);
  feasst_deserialize(vmd_, istr);
}

}  // namespace feasst
