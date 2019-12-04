#include "steppers/include/check_properties.h"

namespace feasst {

class MapCheckProperties {
 public:
  MapCheckProperties() {
    CheckProperties().deserialize_map()["CheckProperties"] = MakeCheckProperties();
  }
};

static MapCheckProperties mapper_check_properties_ = MapCheckProperties();

CheckProperties::CheckProperties(const argtype &args)
  : ModifyUpdateOnly(args) {}

void CheckProperties::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  // store a copy of the system/config
  Configuration config = deep_copy(system->configuration());

  // update properties
  system->unoptimized_energy();

  // see if any eik were not updated properly
  const Select& selection = config.selection_of_all();
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site1 = config.select_particle(part_index).site(site_index);
      const Site& site2 = system->configuration().select_particle(part_index).site(site_index);
      if (!site1.properties().is_equal(site2.properties())) {
        INFO(site1.properties().str());
        INFO(site2.properties().str());
        ERROR("The site properties (old/new) are not equal for " <<
              "part_index:" << part_index << " and site_index:" << site_index);
      }
    }
  }
}

void CheckProperties::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(281, ostr);
}

CheckProperties::CheckProperties(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 281, "version: " << version);
}

}  // namespace feasst
