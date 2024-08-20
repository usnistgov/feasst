#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_inner.h"
#include "system/include/visit_model_cutoff_outer.h"

namespace feasst {

VisitModelCutoffOuter::VisitModelCutoffOuter(argtype * args) : VisitModel() {
//VisitModelCutoffOuter::VisitModelCutoffOuter(argtype * args) : VisitModel(args) {
  // HWH: Strange error if using VisitModel constructor.
  class_name_ = "VisitModelCutoffOuter";
  energy_cutoff_ = dble("energy_cutoff", args, -1);
  if (energy_cutoff_ != -1) {
    ASSERT(energy_cutoff_ > 1e10, "energy_cutoff:" << energy_cutoff_ <<
      " should be > 1e10 to avoid any trial with a chance of being accepted.");
  }
}
VisitModelCutoffOuter::VisitModelCutoffOuter(argtype args) : VisitModelCutoffOuter(&args) {
  feasst_check_all_used(args);
}
VisitModelCutoffOuter::VisitModelCutoffOuter(std::shared_ptr<VisitModelInner> inner,
  argtype args) : VisitModelCutoffOuter(args) {
  set_inner(inner);
}

void VisitModelCutoffOuter::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("visiting model");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  const Select& select_all = config->group_select(group_index);
  bool is_old_config = false;
  if (selection.trial_state() == 0 ||
      selection.trial_state() == 2) {
    is_old_config = true;
  }

  // If possible, query energy map of old configuration instead of pair loop
  if (is_old_config) {
    if (selection.num_particles() == 1) {
      if (get_inner_()->is_energy_map_queryable()) {
        get_inner_()->query_ixn(selection);
        set_energy(inner().energy());
        return;
      }
    }
  }

  // If only one particle in selection, simply exclude part1==part2
  if (selection.num_particles() == 1) {
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      TRACE("part1_index " << part1_index << " s " <<
            selection.particle_indices().size() << " " <<
            selection.site_indices().size());
      for (int select2_index = 0;
           select2_index < select_all.num_particles();
           ++select2_index) {
        const int part2_index = select_all.particle_index(select2_index);
        if (part1_index != part2_index) {
          get_inner_()->set_skip_particle(false);
          for (const int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            for (const int site2_index : select_all.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              if (!inner().skip_particle()) {
                get_inner_()->compute(part1_index, site1_index,
                                      part2_index, site2_index,
                                      config, model_params, model,
                                      is_old_config,
                                      relative_.get(), pbc_.get());
                if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                  set_energy(inner().energy());
                  return;
                }
              }
            }
          }
        }
      }
    }
  } else if (selection.is_equal(config->selection_of_all())) {
    compute_between_selection(model, model_params, selection,
      config, is_old_config, relative_.get(), pbc_.get());

  // If selection is more than one particle but not all particles, skip those in selection
  // Calculate energy in two separate loops.
  } else {
    TRACE("more than one particle in selection");
    for (int select2_index = 0;
         select2_index < select_all.num_particles();
         ++select2_index) {
      const int part2_index = select_all.particle_index(select2_index);
      if (!find_in_list(part2_index, selection.particle_indices())) {
        for (int select1_index = 0;
             select1_index < selection.num_particles();
             ++select1_index) {
          const int part1_index = selection.particle_index(select1_index);
          TRACE("part1_index " << part1_index << " s " <<
                selection.particle_indices().size() << " " <<
                selection.site_indices().size());
          get_inner_()->set_skip_particle(false);
          for (const int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            for (const int site2_index : select_all.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              if (!inner().skip_particle()) {
                get_inner_()->compute(part1_index, site1_index,
                                      part2_index, site2_index,
                                      config, model_params, model,
                                      is_old_config,
                                      relative_.get(), pbc_.get());
                if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                  set_energy(inner().energy());
                  return;
                }
              }
            }
          }
        }
      }
    }

    // In the second loop, compute interactions between different particles in select.
    compute_between_selection(model, model_params, selection,
      config, is_old_config, relative_.get(), pbc_.get());
  }
  set_energy(inner().energy());
}

void VisitModelCutoffOuter::compute_between_selection(
  ModelTwoBody * model,
  const ModelParams& model_params,
  const Select& selection,
  Configuration * config,
  const bool is_old_config,
  Position * relative,
  Position * pbc) {
  for (int select1_index = 0;
       select1_index < selection.num_particles() - 1;
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    TRACE("sel1 " << select1_index << " part1_index " << part1_index << " s "
          << selection.particle_indices().size() << " " <<
          selection.site_indices().size());
    for (int select2_index = select1_index + 1;
         select2_index < selection.num_particles();
         ++select2_index) {
      const int part2_index = selection.particle_index(select2_index);
      if (part1_index != part2_index) {
        TRACE("sel2 " << select2_index << " part2_index " << part2_index);
        get_inner_()->set_skip_particle(false);
        for (const int site1_index : selection.site_indices(select1_index)) {
          TRACE("site1_index " << site1_index);
          for (const int site2_index : selection.site_indices(select2_index)) {
            TRACE("index: " << part1_index << " " << part2_index << " " <<
                  site1_index << " " << site2_index);
            if (!inner().skip_particle()) {
              get_inner_()->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    relative_.get(), pbc_.get());
              if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                set_energy(inner().energy());
                return;
              }
            }
          }
        }
      }
    }
  }
}


class MapVisitModelCutoffOuter {
 public:
  MapVisitModelCutoffOuter() {
    auto obj = MakeVisitModelCutoffOuter();
    obj->deserialize_map()["VisitModelCutoffOuter"] = obj;
  }
};

static MapVisitModelCutoffOuter mapper_ = MapVisitModelCutoffOuter();

VisitModelCutoffOuter::VisitModelCutoffOuter(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2081 == version, version);
  feasst_deserialize(&energy_cutoff_, istr);
}

void VisitModelCutoffOuter::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(2081, ostr);
  feasst_serialize(energy_cutoff_, ostr);
}

}  // namespace feasst
