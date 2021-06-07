#include "utils/include/serialize.h"
#include "chain/include/perturb_reptate.h"

namespace feasst {

class MapPerturbReptate {
 public:
  MapPerturbReptate() {
    auto obj = MakePerturbReptate();
    obj->deserialize_map()["PerturbReptate"] = obj;
  }
};

static MapPerturbReptate mapper_ = MapPerturbReptate();

std::shared_ptr<Perturb> PerturbReptate::create(std::istream& istr) const {
  return std::make_shared<PerturbReptate>(istr);
}

PerturbReptate::PerturbReptate(std::istream& istr)
  : PerturbDistance(istr) {
  ASSERT(class_name_ == "PerturbReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(863 == version, "mismatch version: " << version);
}

void PerturbReptate::serialize_perturb_reptate_(std::ostream& ostr) const {
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(863, ostr);
}

void PerturbReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_reptate_(ostr);
}

void PerturbReptate::finalize(System * system) {
  //PerturbDistance::finalize(system);
  // HWH could also use revert_select instead of finalize_select?
  const Select& mobile = finalize_select()->mobile();
  const int part_index = mobile.particle_indices()[0];
  const Particle& part = system->configuration().select_particle(part_index);
  Select entire(part_index, part);
  Configuration * config = system->get_configuration();
  if (mobile.site_indices()[0][0] == 0) {
    const int site_type = part.site(0).type();
    for (int site = 1; site < entire.num_sites(); ++site) {
      entire.set_site_position(0, site-1, entire.site_positions()[0][site]);
      entire.set_site_properties(0, site-1, entire.site_properties()[0][site]);
      config->set_site_type(part.type(), site - 1, part.site(site).type());
      //config->get_particles_()->get_particle(part_index)->get_site(site-1)->set_cell(
    }
    entire.set_site_position(0, entire.num_sites()-1, part.site(0).position());
    entire.set_site_properties(0, entire.num_sites()-1, part.site(0).properties());
    config->set_site_type(part.type(), entire.num_sites() - 1, site_type);
  } else {
    const int site_type = part.site(entire.num_sites()-1).type();
    for (int site = entire.num_sites() - 1; site >= 1; --site) {
      entire.set_site_position(0, site, entire.site_positions()[0][site - 1]);
      entire.set_site_properties(0, site, entire.site_properties()[0][site - 1]);
      config->set_site_type(part.type(), site, part.site(site - 1).type());
    }
    entire.set_site_position(0, 0, part.site(part.num_sites()-1).position());
    entire.set_site_properties(0, 0, part.site(part.num_sites()-1).properties());
    config->set_site_type(part.type(), 0, site_type);
  }
  // HWH this may no longer be necessary now that cells aren't in config
  config->update_positions(entire, false);
}

}  // namespace feasst
