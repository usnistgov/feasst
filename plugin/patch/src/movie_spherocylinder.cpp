#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "patch/include/movie_spherocylinder.h"

namespace feasst {

FEASST_MAPPER(MovieSpherocylinder, argtype({{"output_file", "place_holder"}}));

MovieSpherocylinder::MovieSpherocylinder(argtype * args) : AnalyzeWriteOnly(args) {
  set_append();
  ASSERT(!output_file().empty(), "file name is required");
  args->insert({"append", "true"}); // always append
  xyz_ = FileXYZSpherocylinder(args);
}
MovieSpherocylinder::MovieSpherocylinder(argtype args) : MovieSpherocylinder(&args) { feasst_check_all_used(args); }

void MovieSpherocylinder::initialize(MonteCarlo * mc) {
  const Criteria& criteria = mc->criteria();
  const std::string name = output_file(criteria);
  ASSERT(!name.empty(), "file name required. Did you forget to " <<
    "Analyze::set_output_file()?");

  // write xyz
  if (state() == criteria.state()) {
    xyz_.write(name, configuration(mc->system()));
  }

  // write vmd
  std::stringstream ss;
  ss << name << ".vmd";
  vmd_.write(ss.str(), configuration(mc->system()), name);
}

std::string MovieSpherocylinder::write(const MonteCarlo& mc) {
  // ensure the following order matches the header from initialization.
  xyz_.write(output_file(mc.criteria()), configuration(mc.system()));
  return std::string("");
}

void MovieSpherocylinder::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2364, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
  feasst_serialize_fstobj(vmd_, ostr);
}

MovieSpherocylinder::MovieSpherocylinder(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2364, "version mismatch:" << version);
  feasst_deserialize_fstobj(&xyz_, istr);
  feasst_deserialize_fstobj(&vmd_, istr);
}

}  // namespace feasst
