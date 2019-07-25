#include <vector>
#include "system/include/visit_model_intra.h"
#include "system/include/model_two_body.h"

namespace feasst {

void VisitModelIntra::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  TRACE("intra particle energy_of_selection");
  ASSERT(group_index == 0, "need to implement site1 loop filtering particles by group");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_);
  for (int sp1index = 0;
       sp1index < static_cast<int>(selection.particle_indices().size());
       ++sp1index) {
    const int part1_index = selection.particle_index(sp1index);
  //for (int part1_index : selection.particle_indices()) {
    TRACE("particle: " << part1_index);
    const Particle& part1 = config->select_particle(part1_index);
    // the first site loop is over all sites in part1 and group_index
    // the second is all sites in selection
    // HWH optimize this

    // here we use excluded to account for chain regrowth, etc.
    // exclude the particles which haven't been grown yet.
    // or exclude particles which will form new bonds (reptate).
    Select sites1 = selection;
    if (selection.excluded()) {
      sites1.remove(*(selection.excluded()));
      TRACE("excluded " << selection.excluded()->str());
    }
    TRACE("sites1: " << sites1.str());
    const std::vector<int>& site1_indices = sites1.site_indices(0);

    Select sites2;
    sites2.add_particle(part1, part1_index);
    if (selection.excluded()) {
      sites2.remove(*(selection.excluded()));
      TRACE("excluded " << selection.excluded()->str());
    }
    TRACE("sites2: " << sites2.str());
    const std::vector<int>& site2_indices = sites2.site_indices(0);

    for (int sel1_index = 0;
         sel1_index < static_cast<int>(site1_indices.size());
         ++sel1_index) {
      TRACE("sel1_index " << sel1_index << " size " << site1_indices.size());
      const int site1_index = site1_indices[sel1_index];
      for (int sel2_index = 0;
           sel2_index < static_cast<int>(site2_indices.size());
           ++sel2_index) {
        const int site2_index = site2_indices[sel2_index];

        // if sites in particle selection > 1, attempt the following check.
        // if site2 is in selection, then require site1 < site2
        if (!find_in_list(site2_index, site1_indices) or
            site1_index < site2_index) {

          // here we determine if the pair of sites is forced to be included
          bool include = false;
          if (selection.old_bond()) {
            const int incl1_site = site1_indices[0];
            const int incl2_site = selection.old_bond()->site_indices()[sp1index][0];
            if ( (site1_index == incl1_site and
                  site2_index == incl2_site) or
                 (site1_index == incl2_site and
                  site2_index == incl1_site) ) {
              include = true;
              TRACE("include " << include << " incl " << incl1_site << " " << incl2_site);
            }
          }

          // here we determine if the pair of sites is forced to be excluded
          bool exclude = false;
          if (selection.new_bond()) {
            const int excl1_site = site1_indices[0];
            const int excl2_site = selection.new_bond()->site_indices()[sp1index][0];
            if ( (site1_index == excl1_site and
                  site2_index == excl2_site) or
                 (site1_index == excl2_site and
                  site2_index == excl1_site) ) {
              exclude = true;
              TRACE("exclude " << exclude << " excl " << excl1_site << " " << excl2_site);
            }
          }

          // forced exclude takes precedent over forced include
          if ( (include or std::abs(site1_index - site2_index) > intra_cut_) and (!exclude) ) {
            TRACE("sites: " << site1_index << " " << site2_index);
            inner()->compute(part1_index, site1_index, part1_index, site2_index,
                             config, model_params, model, &relative_);
          }
        }
      }
    }
  }
  set_energy(inner()->energy());
}

class MapVisitModelIntra {
 public:
  MapVisitModelIntra() {
    VisitModelIntra().deserialize_map()["VisitModelIntra"] =
      std::make_shared<VisitModelIntra>();
  }
};

static MapVisitModelIntra mapper_ = MapVisitModelIntra();

std::shared_ptr<VisitModel> VisitModelIntra::create(std::istream& istr) const {
  return std::make_shared<VisitModelIntra>(istr);
}

VisitModelIntra::VisitModelIntra(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(754 == version, version);
  feasst_deserialize(&intra_cut_, istr);
}

void VisitModelIntra::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(754, ostr);
  feasst_serialize(intra_cut_, ostr);
}



}  // namespace feasst
