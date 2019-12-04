#include "steppers/include/movie.h"

namespace feasst {

class MapMovie {
 public:
  MapMovie() {
    auto obj = MakeMovie({{"file_name", "place_holder"}});
    obj->deserialize_map()["Movie"] = obj;
  }
};

static MapMovie mapper_ = MapMovie();

Movie::Movie(const argtype &args ) : AnalyzeWriteOnly(args) {
  set_append();

  // require file_name
  args_.init(args);
  ASSERT(args_.key("file_name").used(), "file name is required");
}

void Movie::initialize(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ASSERT(!file_name().empty(), "file name required. Did you forget to " <<
    "Analyze::set_file_name()?");

  // write xyz
  xyz_.set_append(1);
  if (state() == criteria->state()) {
    xyz_.write(file_name(), system.configuration());
  }

  // write vmd
  std::stringstream ss;
  ss << file_name() << ".vmd";
  vmd_.write(ss.str(), system.configuration(), file_name());
}

std::string Movie::write(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  xyz_.write(file_name(), system.configuration());
  return std::string("");
}

void Movie::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(536, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
  feasst_serialize_fstobj(vmd_, ostr);
}

Movie::Movie(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 536, "version mismatch:" << version);
  feasst_deserialize_fstobj(&xyz_, istr);
  feasst_deserialize_fstobj(&vmd_, istr);
}

}  // namespace feasst
