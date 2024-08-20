#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "netcdf/include/movie_netcdf.h"

namespace feasst {

class MapMovieNETCDF {
 public:
  MapMovieNETCDF() {
    auto obj = MakeMovieNETCDF();
    //auto obj = MakeMovieNETCDF({{"file_name", feasst::install_dir() + "/plugin/netcdf/test/data/first.netCDF.nc"}});
    //auto obj = MakeMovieNETCDF({{"file_name", "place_holder"}});
    obj->deserialize_map()["MovieNETCDF"] = obj;
  }
};

static MapMovieNETCDF mapper_ = MapMovieNETCDF();

MovieNETCDF::MovieNETCDF(argtype * args)
: AnalyzeWriteOnly(args),
  xyz_(file_parse_(args)) {
  set_append();
}
MovieNETCDF::MovieNETCDF(argtype args) : MovieNETCDF(&args) {
  feasst_check_all_used(args);
}

void MovieNETCDF::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const std::string name = file_name(*criteria);
  ASSERT(!name.empty(), "file name required. Did you forget to " <<
    "Analyze::set_file_name()?");
  xyz_->initialize(system->configuration());
  // write xyz
  if (state() == criteria->state()) {
    xyz_->write(system->configuration());
  }
}

std::string MovieNETCDF::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  xyz_->write(system.configuration());
  //xyz_.write(file_name(criteria), system.configuration());
  return std::string("");
}

void MovieNETCDF::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7698, ostr);
  //feasst_serialize_fstobj(xyz_, ostr);
}

MovieNETCDF::MovieNETCDF(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7698, "version mismatch:" << version);
  //feasst_deserialize_fstobj(&xyz_, istr);
}

}  // namespace feasst
