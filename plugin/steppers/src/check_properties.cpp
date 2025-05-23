#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/check_properties.h"

namespace feasst {

FEASST_MAPPER(CheckProperties,);

CheckProperties::CheckProperties(argtype * args) : ModifyUpdateOnly(args) {
  tolerance_ = dble("tolerance", args, 1e-15);
}
CheckProperties::CheckProperties(argtype args) : CheckProperties(&args) {
  feasst_check_all_used(args);
}

void CheckProperties::update(MonteCarlo * mc) {
  System * system = mc->get_system();
  for (int iconf = 0; iconf < system->num_configurations(); ++iconf) {
    // store a copy of the system/config
    Configuration config = deep_copy(system->configuration(iconf));

    // update properties
    system->unoptimized_energy(iconf);

    // see if any properties do not match the copy
    bool failure = false;
    const Select& selection = config.selection_of_all();
    for (int select_index = 0;
         select_index < selection.num_particles();
         ++select_index) {
      const int part_index = selection.particle_index(select_index);
      for (int site_index : selection.site_indices(select_index)) {
        const Site& site1 = config.select_particle(part_index).site(site_index);
        const Site& site2 = system->configuration(iconf).select_particle(part_index).site(site_index);
        if (!site1.properties().is_equal(site2.properties(), tolerance_)) {
          INFO(site1.properties().str());
          INFO(site2.properties().str());
          INFO("The site properties (old/new) are not equal for " <<
               "part_index:" << part_index << " and site_index:" << site_index);
          failure = true;
        }
      }
      ASSERT(!failure, "check failed");
    }
  }
}

void CheckProperties::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(281, ostr);
  feasst_serialize(tolerance_, ostr);
}

CheckProperties::CheckProperties(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 281, "version: " << version);
  feasst_deserialize(&tolerance_, istr);
}

}  // namespace feasst
