#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/configuration.h"
#include "system/include/bond_visitor.h"

namespace feasst {

class MapBondVisitor {
 public:
  MapBondVisitor() {
    auto obj = MakeBondVisitor();
    obj->deserialize_map()["BondVisitor"] = obj;
  }
};

static MapBondVisitor mapper_ = MapBondVisitor();

std::map<std::string, std::shared_ptr<BondVisitor> >& BondVisitor::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondVisitor> >* ans =
     new std::map<std::string, std::shared_ptr<BondVisitor> >();
  return *ans;
}

BondVisitor::BondVisitor(argtype args) {
  verbose_ = boolean("verbose", &args, false);
  if (VERBOSE_LEVEL == 5) verbose_ = true;
  FEASST_CHECK_ALL_USED(args);
}

void BondVisitor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_bond_visitor_(ostr);
}

std::shared_ptr<BondVisitor> BondVisitor::create(std::istream& istr) const {
  return std::make_shared<BondVisitor>(istr);
}

std::shared_ptr<BondVisitor> BondVisitor::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondVisitor::serialize_bond_visitor_(std::ostream& ostr) const {
  feasst_serialize_version(303, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize(verbose_, ostr);
}

BondVisitor::BondVisitor(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(303 == version, "mismatch version: " << version);
  feasst_deserialize(&energy_, istr);
  feasst_deserialize(&verbose_, istr);
}

void BondVisitor::compute_two(const Select& selection,
    const Configuration& config) {
  double en = 0.;
  //TRACE(selection.str());
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const Particle& unique_part = config.unique_type(part.type());
    const Particle& part_type = config.particle_type(part.type());
    if (part_type.num_bonds() == 0) continue;
    for (int site0_index : selection.site_indices(select_index)) {
      DEBUG("site0_index " << site0_index);
      const Site& site0 = part.site(site0_index);
      for (int site1_index : part_type.bond_neighbors(site0_index)) {
        DEBUG("site1_index " << site1_index);
        const Site& site1 = part.site(site1_index);
        if (site1.is_physical()) {
          if (site0_index < site1_index ||
              !find_in_list(site1_index,
                            selection.site_indices(select_index))) {
            const Position& ri = site0.position();
            const Position& rj = site1.position();
            const Bond& bond_type = part_type.bond(site0_index, site1_index);
            const Bond& bond = unique_part.bond(bond_type.type());
            ASSERT(bond_.deserialize_map().count(bond.model()) == 1,
              "bond model " << bond.model() << " not recognized.");
            en += bond_.deserialize_map()[bond.model()]->energy(
              ri, rj, bond);
            if (verbose_) {
              if (std::abs(en) > NEAR_ZERO) {
                TRACE("bond ij " << part_index << " " << bond.site(0) << " "
                  << bond.site(1) << " sq " << ri.squared_distance(rj));
              }
            }
          }
        }
      }
    }
  }
  energy_two_body_ = en;
}

void BondVisitor::compute_three(
    const Select& selection,
    const Configuration& config) {
  TRACE("BondThreeBody of " << selection.str());
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const Particle& part_type = config.particle_type(part.type());
    if (part_type.num_angles() == 0) continue;
    const Particle& unique_part = config.unique_type(part.type());
    for (int site0_index : selection.site_indices(select_index)) {
      const Site& site0 = part.site(site0_index);
      DEBUG("site0_index " << site0_index);
      for (const std::vector<int>& ang : part_type.angle_neighbors(site0_index)) {
        const int site1_index = ang[0];
        const int site2_index = ang[1];
        const Site& site1 = part.site(site1_index);
        const Site& site2 = part.site(site2_index);
        if (site1.is_physical() && site2.is_physical()) {
          /*
             three bodies are defined by sites 0-1-2, which also appear as 2-1-0
             but never 1-0-2, as that would change the angle/vertex.
             Therefore, when avoiding double counts for whole-molecule compute,
             check that site0 < site2, unless site2 isn't in selection.
           */
          if (site0_index < site2_index ||
              !find_in_list(site2_index, selection.site_indices(select_index))) {
            TRACE("sites " << site0_index << " " << site1_index << " " << site2_index);
            const Position * ri = &site0.position();
            const Position * rj = &site1.position();
            const Position * rk = &site2.position();
            const Angle& angle_type = part_type.angle(site0_index,
              site1_index, site2_index);
            const Angle& angle = unique_part.angle(angle_type.type());
            ASSERT(angle_.deserialize_map().count(angle.model()) == 1,
              "angle model " << angle.model() << " not recognized.");

            // In 2D, angle i-j-k is not the same as k-j-i.
            // angle i-j-k = 2pi - angle k-j-i
            if (ri->dimension() == 2) {
              TRACE("site0 of angle " << angle.site(0));
              if (angle.site(0) != site0_index) {
                ri = &site2.position();
                rk = &site0.position();
              }
            }
            en += angle_.deserialize_map()[angle.model()]->energy(
              *ri, *rj, *rk, angle);
            if (verbose_) {
              if (std::abs(en) > NEAR_ZERO) {
                const double ang = rj->vertex_angle_radians(*ri, *rk);
                const double theta0 = degrees_to_radians(angle.property("theta0"));
                TRACE("angle " << part_index << " ijk " << angle.site(0) << " "
                  << angle.site(1) << " " << angle.site(2) << " "
                  << "ri " << ri->str() << " "
                  << "rj " << rj->str() << " "
                  << "rk " << rk->str() << " "
                  << "ang " << ang << " "
                  << "theta0 " << theta0 << " "
                  << "diff " << ang - theta0);
              }
            }
          }
        }
      }
    }
  }
  energy_three_body_ = en;
}

void BondVisitor::compute_four(
    const Select& selection,
    const Configuration& config) {
  double en = 0.;
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(part_index);
    const Particle& unique_part = config.unique_type(part.type());
    const Particle& part_type = config.particle_type(part.type());
    if (part_type.num_dihedrals() == 0) continue;
    for (int site0_index : selection.site_indices(select_index)) {
      const Site& site0 = part.site(site0_index);
      for (const std::vector<int>& dih : part_type.dihedral_neighbors(site0_index)) {
        const int site1_index = dih[0];
        const int site2_index = dih[1];
        const int site3_index = dih[2];
        const Site& site1 = part.site(site1_index);
        const Site& site2 = part.site(site2_index);
        const Site& site3 = part.site(site3_index);
        if (site1.is_physical() && site2.is_physical() && site3.is_physical()) {
          /*
             four bodies defined by sites 0-1-2-3 also appear as 3-2-1-0,
             but never 0-2-1-3, etc, as that would change the angle.
             Therefore, when avoiding double counts for whole-molecule compute,
             check that site0 < site3, unless site3 isn't in selection.
           */
          if (site0_index < site3_index ||
              !find_in_list(site3_index, selection.site_indices(select_index))) {
            const Position& ri = site0.position();
            const Position& rj = site1.position();
            const Position& rk = site2.position();
            const Position& rl = site3.position();
            TRACE("ri " << ri.str());
            TRACE("rj " << rj.str());
            TRACE("rk " << rk.str());
            TRACE("rl " << rl.str());
            const Dihedral& dihedral_type = part_type.dihedral(site0_index, site1_index,
              site2_index, site3_index);
            TRACE("type of dihedral " << dihedral_type.type());
            const Dihedral& dihedral = unique_part.dihedral(dihedral_type.type());
            TRACE("model of dihedral " << dihedral.model());
            ASSERT(dihedral_.deserialize_map().count(dihedral.model()) == 1,
              "dihedral model " << dihedral.model() << " not recognized.");
            en += dihedral_.deserialize_map()[dihedral.model()]->energy(
              ri, rj, rk, rl, dihedral);
            TRACE("en: " << en << " sites " << site0_index << " " << site1_index << " " << site2_index << " " << site3_index);
            if (verbose_) {
              if (std::abs(en) > NEAR_ZERO) {
                FATAL("not impl");
              }
            }
          }
        }
      }
    }
  }
  energy_four_body_ = en;
}

void BondVisitor::compute_all(const Select& selection,
                              const Configuration& config) {
  compute_two(selection, config);
  compute_three(selection, config);
  compute_four(selection, config);
  energy_ = energy_two_body_ + energy_three_body_ + energy_four_body_;
}

void BondVisitor::compute_all(const Configuration& config,
                              const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_all(selection, config);
}

void BondVisitor::compute_two(const Configuration& config,
                              const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_two(selection, config);
}

void BondVisitor::compute_three(const Configuration& config,
                                const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_three(selection, config);
}

void BondVisitor::compute_four(const Configuration& config,
                               const int group_index) {
  const Select& selection = config.group_selects()[group_index];
  compute_four(selection, config);
}

}  // namespace feasst
