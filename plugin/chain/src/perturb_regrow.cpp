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
  DEBUG("bonded_to " << site_bonded_to << " to update " << site_to_update);
  const Bond& bond = selection().bond(0, site_to_update, site_bonded_to, config);
  const double l0 = bond.property("l0");
//  DEBUG(
//    "to update " << site_to_update << " "
//    "bonded to " << site_bonded_to << " "
//    "l0 " << l0
//  );

  // obtain a new position for the site to update
//  DEBUG("bonded sites " << bonded->num_sites());
  set_in_sphere(l0, index, bonded);
}

void PerturbRegrow::set_in_sphere(const double bond_length, const int index, SelectList * bonded) {
  DEBUG("num sites " << anchor().num_sites());
  const Position& anchor_spot = anchor().site_positions()[0][index];
  // const Position& anchor_spot = anchor().site_positions()[0][index];
  DEBUG("dim " << anchor_spot.dimension());
  Position rebond = anchor_spot;
  random_.unit_sphere_surface(&rebond);
  rebond.multiply(bond_length);
  rebond.add(anchor_spot);
  DEBUG("new position " << rebond.str());
  bonded->set_site_position(0, index, rebond);
}

}  // namespace feasst
