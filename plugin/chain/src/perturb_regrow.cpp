#include "chain/include/perturb_regrow.h"

namespace feasst {

void PerturbRegrow::rebond_(
    const int index,
    System * system,
    SelectList * bonded) {
  // obtain the bond length
  const Configuration& config = system->configuration();
  const int site_bonded_to = anchor().site_indices(0)[index];
  const int site_to_update = selection().site_indices(0)[index];
  const Bond& bond = selection().bond(0, site_to_update, site_bonded_to, config);
  const double l0 = bond.property("l0");
//  INFO(
//    "to update " << site_to_update << " "
//    "bonded to " << site_bonded_to << " "
//    "l0 " << l0
//  );
  // obtain a new position for the site to update
//  INFO("bonded sites " << bonded->num_sites());
  const Position& anchor_spot = anchor().site_positions()[0][index];
  // const Position& anchor_spot = anchor().site_positions()[0][index];
//  INFO("dim " << anchor_spot.dimension());
  Position rebond = anchor_spot;
  random_.unit_sphere_surface(&rebond);
  DEBUG("bonded_to " << site_bonded_to << " to update " << site_to_update);
  rebond.multiply(l0);
  rebond.add(anchor_spot);

  bonded->set_site_position(0, index, rebond);
}

}  // namespace feasst
