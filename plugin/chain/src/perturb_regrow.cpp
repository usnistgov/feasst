#include "chain/include/perturb_regrow.h"

namespace feasst {

void PerturbRegrow::rebond_(
    const int site_to_update,
    const int site_bonded_to,
    System * system,
    SelectList * bonded) {
  // obtain the bond length
  const Configuration& config = system->configuration();
  const Bond& bond = selection().bond(0, site_to_update, site_bonded_to, config);
  const double l0 = bond.property("l0");

  // obtain a new position for the site to update
  const Position& anchor = bonded->site_positions()[0][site_bonded_to];
  Position rebond = anchor;
  random_.unit_sphere_surface(&rebond);
  rebond.multiply(l0);
  rebond.add(anchor);

  bonded->set_site_position(0, site_to_update, rebond);
}

}  // namespace feasst
