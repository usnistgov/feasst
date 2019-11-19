
#ifndef FEASST_STEPPERS_CHECK_PROPERTIES_H_
#define FEASST_STEPPERS_CHECK_PROPERTIES_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Check equivalence of stored properties and updated properties.
 */
class CheckProperties : public ModifyUpdateOnly {
 public:
  CheckProperties(const argtype &args = argtype()) : ModifyUpdateOnly(args) {}

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override {
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

  std::string class_name() const override { return std::string("CheckProperties"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(281, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CheckProperties>(istr); }

  CheckProperties(std::istream& istr) : ModifyUpdateOnly(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(version == 281, "version: " << version);
  }

 private:
  double tolerance_;
};

inline std::shared_ptr<CheckProperties> MakeCheckProperties(const argtype &args = argtype()) {
  return std::make_shared<CheckProperties>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_PROPERTIES_H_
