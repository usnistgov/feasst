
#ifndef FEASST_CHAIN_PERTURB_H_
#define FEASST_CHAIN_PERTURB_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

// HWH move to chain plugin
// HWH recenter particle position
class PerturbPivot : public PerturbRotate {
 public:
  PerturbPivot(const argtype& args = argtype()) : PerturbRotate(args) {}

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot to the anchor.
  /// Dont rotate the particle positions.
  void move(System * system,
      TrialSelect * select) override {
    const Position& pivot = select->anchor_position(0, 0, system);
    DEBUG("piv " << pivot.str());
    PerturbRotate::move(system, select, pivot, false);
    DEBUG(select->mobile().site_positions()[0][0].str());
    DEBUG(select->mobile().site_positions()[0][1].str());
  }
};

// HWH move to chain plugin
// HWH recenter particle position
class PerturbCrankshaft : public PerturbRotate {
 public:
  PerturbCrankshaft(const argtype& args = argtype()) : PerturbRotate(args) {}

  /// Set the pivot and axis of rotation by the ends of the selection.
  /// Select rotation angle randomly, bounded by tunable parameter.
  /// Dont rotate the particle positions.
  void move(System * system,
      TrialSelect * select) override {
    const Position& pivot = select->mobile().site_positions()[0].front();
    axis_ = select->mobile().site_positions()[0].back();
    axis_.subtract(pivot);
    axis_.normalize();
    const double max_angle = tunable().value();
    const double angle = random()->uniform_real(-max_angle, max_angle);
    rot_mat_.axis_angle(axis_, angle);
    PerturbRotate::move(pivot, rot_mat_, system, select,
      false // do not rotate particle positions
    );
  }

 private:
  // temporary
  Position axis_;
  RotationMatrix rot_mat_;
};

/**
  For a reptation, if new bond is accepted, then change the positions of all the
  sites along the chain.
 */
class PerturbReptate : public PerturbDistanceFromAnchor {
 public:
  PerturbReptate(const argtype& args = argtype()) : PerturbDistanceFromAnchor(args) {}
  void move(System * system,
      TrialSelect * select) override {
    PerturbDistanceFromAnchor::move(system, select);
    set_finalize_possible(true, select);
  }

  void finalize(System * system) override {
    const SelectList& mobile = finalize_select()->mobile();
    const int part_index = mobile.particle_indices()[0];
    SelectList entire = SelectList().particle(part_index,
                                              system->configuration(),
                                              0 // group that includes all
                                              );
    if (mobile.site_indices()[0][0] == 0) {
      for (int site = 1; site < entire.num_sites(); ++site) {
        entire.set_site_position(0, site - 1, entire.site_positions()[0][site]);
        entire.set_site_properties(0, site - 1, entire.site_properties()[0][site]);
      }
      entire.set_site_position(0, entire.num_sites() - 1, mobile.site_positions()[0][0]);
      entire.set_site_properties(0, entire.num_sites() - 1, mobile.site_properties()[0][0]);
    } else {
      for (int site = entire.num_sites() - 1; site >= 1; --site) {
        entire.set_site_position(0, site, entire.site_positions()[0][site - 1]);
        entire.set_site_properties(0, site, entire.site_properties()[0][site - 1]);
      }
      entire.set_site_position(0, 0, mobile.site_positions()[0][0]);
      entire.set_site_properties(0, 0, mobile.site_properties()[0][0]);
    }
    DEBUG("entire " << entire.str() << " pos " << entire.site_positions()[0][0].str() << " end " << entire.site_positions()[0][49].str());
    system->get_configuration()->update_positions(entire, false);
  }
};

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_H_
