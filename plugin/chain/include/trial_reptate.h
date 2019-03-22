
#ifndef FEASST_CHAIN_TRIAL_REPTATE_H_
#define FEASST_CHAIN_TRIAL_REPTATE_H_

#include "core/include/trial.h"
#include "core/include/utils_io.h"
#include "chain/include/perturb_regrow.h"

namespace feasst {

/**
  HWH bug where the intracut ignores interaction of new end with the particle
  bonded to the old end.
  For example we may need site 0 to interact with site 1 but bond connectivity
  failes

  Proposed fixes:
  4. update positions of all, but compute energy only for selection

  3. change the site_index in selection ? ... you'd get the wrong positions

  2. update positions of all sites in the chain (slowwww)

  1. introduce new particle types 1-2-...-49-0 bond connectivity and
      49-0-1-...-47-48 bond connectivity.
      temporarily change the type while computing energies of the selection
      if accepted, revert type and update entire chain
*/

/**
  Select random particle
  take one of the end points and attempt to attach it to the other end point.
  This would involve using a PerturbRegrow to attach the endpoint.

  Compute the energy of the old end point.
  Minus the energy of the new end point.

  But the selection contains the entire particle which needs to be updated.
 */
class TrialReptate : public Trial {
 public:
  TrialReptate(const argtype& args = argtype()) : Trial(args) {}

  /// Refactor to attempt, including perturb-revert or acceptance with new selection
  void attempt(Criteria * criteria, System * system) override {
    before_attempt(criteria, system, &regrow_, &accept_criteria_);
    Configuration * config = regrow_.get_config_before_move(system);
    const Configuration& const_config = system->configuration();
    const int max_length = 1;
    SelectList select;
    select.random_end_segment_in_particle(group_index(), const_config, max_length);
    regrow_.set_selection(select);
    // const double bond_length = select_and_bond_length(system);
    SelectList bonded;
    if (regrow_.selection().is_empty()) {
      accept_criteria_.force_rejection = 1;
    } else {
      const double pe_old = system->energy(regrow_.selection());

      // set the anchor as the other end segment
      SelectList anchor;
      const int part_index = regrow_.selection().particle_indices()[0];
      int site_bonded_to = -1;
      const int site_to_update = regrow_.selection().site_indices()[0][0];
      int other_bond = -1;
      const int particle_index = regrow_.selection().particle_indices()[0];
      const int num_sites = const_config.select_particle(particle_index).num_sites();
      if (site_to_update == 0) {
        site_bonded_to = num_sites - 1;
        other_bond = 1;
      } else {
        site_bonded_to = 0;
        other_bond = num_sites - 2;
      }
      anchor.add_site(part_index, site_bonded_to);
      anchor.resize();
      anchor.load_positions(const_config.particles());
      regrow_.set_anchor(anchor);

      // exclude the anchor
      select.set_new_bond(anchor);

      // include the the particle that use to be bonded to select
      SelectList previous;
      previous.add_site(part_index, other_bond);
      previous.resize();
      previous.load_positions(const_config.particles());
      select.set_old_bond(previous);

      // set the selection
      regrow_.set_selection(select);
      regrow_.set_selection_state("old");

      DEBUG("site_to_update " << site_to_update << " site_bonded_to " << site_bonded_to << " other_bond " << other_bond);
      const Bond& bond = anchor.bond(0, site_to_update, other_bond, const_config);
      const double bond_length = bond.property("l0");

      bonded = regrow_.selection();
      regrow_.set_in_sphere(bond_length, 0, &bonded);
      config->update_positions(bonded);
      DEBUG("bonded " << bonded.str() << " pos " << bonded.site_positions()[0][0].str());
      DEBUG("sel " << regrow_.selection().str());
      regrow_.after_move();

      const double pe_new = system->energy(regrow_.selection());
      DEBUG("pe_new " << pe_new);
      const double delta_energy = pe_new - pe_old;
      accept_criteria_.ln_metropolis_prob += -criteria->beta()*delta_energy;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.energy_new_select = pe_new;
      accept_criteria_.system = system;
      DEBUG("delta_energy " << delta_energy);
    }

    if (criteria->is_accepted(accept_criteria_)) {
      DEBUG("accepted");
      record_success();

      // update entire particle upon success.
      DEBUG("bonded " << bonded.str() << " pos " << bonded.site_positions()[0][0].str());
      const int part_index = bonded.particle_indices()[0];
      SelectList entire = SelectList().particle(part_index,
                                                system->configuration(),
                                                0 // group that includes all
                                                );
      if (regrow_.selection().site_indices()[0][0] == 0) {
        for (int site = 1; site < entire.num_sites(); ++site) {
          entire.set_site_position(0, site - 1, entire.site_positions()[0][site]);
          entire.set_site_properties(0, site - 1, entire.site_properties()[0][site]);
        }
        entire.set_site_position(0, entire.num_sites() - 1, bonded.site_positions()[0][0]);
        entire.set_site_properties(0, entire.num_sites() - 1, bonded.site_properties()[0][0]);
      } else {
        for (int site = entire.num_sites() - 1; site >= 1; --site) {
          entire.set_site_position(0, site, entire.site_positions()[0][site - 1]);
          entire.set_site_properties(0, site, entire.site_properties()[0][site - 1]);
        }
        entire.set_site_position(0, 0, bonded.site_positions()[0][0]);
        entire.set_site_properties(0, 0, bonded.site_properties()[0][0]);
      }
      DEBUG("entire " << entire.str() << " pos " << entire.site_positions()[0][0].str() << " end " << entire.site_positions()[0][49].str());
      config->update_positions(entire);
    } else {
      DEBUG("rejected");
      DEBUG("sel " << regrow_.selection().str());
      DEBUG("pos " << regrow_.selection().site_positions()[0][0].str());
      regrow_.revert();
    }
  }

  virtual ~TrialReptate() {}

 private:
  AcceptanceCriteria accept_criteria_;
  PerturbRegrow regrow_;
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_REPTATE_H_
