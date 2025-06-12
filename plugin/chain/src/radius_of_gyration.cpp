#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/histogram.h"
#include "math/include/position.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "chain/include/radius_of_gyration.h"

namespace feasst {

FEASST_MAPPER(RadiusOfGyration,);

RadiusOfGyration::RadiusOfGyration(argtype * args) : Analyze(args) {
  group_index_ = integer("group_index", args, 0);
  if (boolean("print_histogram", args, false)) {
    hist_ = std::make_shared<Histogram>(args);
  }
}
RadiusOfGyration::RadiusOfGyration(argtype args) : RadiusOfGyration(&args) {
  feasst_check_all_used(args);
}

void RadiusOfGyration::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string RadiusOfGyration::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator().status_header() << ",rgu,rguu" << std::endl;
  return ss.str();
}

void RadiusOfGyration::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  const Select& selection = system.configuration().group_select(group_index_);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = system.configuration().select_particle(part_index);
    Position r_cm(system.configuration().dimension());
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      r_cm.add(site.position());
    }
    r_cm.divide(selection.num_sites());
    double rg = 0.;
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      rg += site.position().squared_distance(r_cm);
    }
    const double rgn = std::sqrt(rg/selection.num_sites());
    get_accumulator()->accumulate(rgn);
    if (hist_) {
      hist_->add(rgn);
    }
    const double en = mc.criteria().current_energy();
    rg_e_.accumulate(rgn*en);
    rg_e2_.accumulate(rgn*en*en);
  }
}

std::string RadiusOfGyration::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << accumulator().status() << "," << rg_e_.average();
  ss << "," << rg_e2_.average();
  ss << std::endl;
  DEBUG(ss.str());
  if (hist_) {
    ss << hist_->str();
    ss << std::endl;
  }
  return ss.str();
}

void RadiusOfGyration::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7685, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstobj(rg_e_, ostr);
  feasst_serialize_fstobj(rg_e2_, ostr);
  feasst_serialize(hist_, ostr);
}

RadiusOfGyration::RadiusOfGyration(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7685, "mismatch version:" << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize_fstobj(&rg_e_, istr);
  feasst_deserialize_fstobj(&rg_e2_, istr);
//  HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(hist_, istr);
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      hist_ = std::make_shared<Histogram>(istr);
    }
  }
}

RadiusOfGyration::RadiusOfGyration(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = RadiusOfGyration(ss);
}

}  // namespace feasst
