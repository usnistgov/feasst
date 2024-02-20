#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "chain/include/perturb_connector.h"

namespace feasst {

PerturbConnector::PerturbConnector(argtype args) : PerturbConnector(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbConnector::PerturbConnector(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbConnector";
}

class MapPerturbConnector {
 public:
  MapPerturbConnector() {
    auto obj = MakePerturbConnector();
    obj->deserialize_map()["PerturbConnector"] = obj;
  }
};

static MapPerturbConnector mapper_ = MapPerturbConnector();

std::shared_ptr<Perturb> PerturbConnector::create(std::istream& istr) const {
  return std::make_shared<PerturbConnector>(istr);
}

void PerturbConnector::move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) {
  if (is_position_held) return;
  ASSERT(select->mobile().num_sites() == 1,
    "selection num_sites: " << select->mobile().num_sites());
  ASSERT(select->mobile().site_positions().size() == 1, "requires coordinates");
  DEBUG("num particle positions " << select->mobile().site_positions().size());
  DEBUG("num site positions " << select->mobile().site_positions()[0].size());
  ASSERT(select->anchor().num_sites() == 1, "selection error");
  const Position& anchor = select->anchor_position(0, 0, *system);
  DEBUG("anchor " << anchor.str());
  Position * mobile = select->get_mobile()->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("dist " << mobile->distance(anchor));

  // set mobile to position of selected site in rigid definition in fstprt.
  const Configuration& config = select->configuration(*system);
  const Particle& part = config.select_particle(select->mobile().particle_index(0));
  const int particle_type = part.type();
  DEBUG("particle_type " << particle_type);
  const int site_type = part.site(select->mobile().site_index(0, 0)).type();
  DEBUG("site_type " << site_type);
  const Particle& part_type = config.particle_type(part.type());
  const Position& fixed_mobile = part_type.site(select->mobile().site_index(0, 0)).position();
  DEBUG("fixed mobile " << fixed_mobile.str());
  const Position& fixed_anchor = part_type.site(select->anchor().site_index(0, 0)).position();
  DEBUG("fixed anchor " << fixed_anchor.str());
  *mobile = fixed_mobile;
  mobile->subtract(fixed_anchor);
  DEBUG("mobile " << mobile->str());

  // rotate mobile by current eulers in anchor
  const Euler& euler = config.select_particle(select->anchor().particle_index(0)).site(select->anchor().site_index(0, 0)).euler();
  DEBUG("euler " << euler.str());
  euler.compute_rotation_matrix(&rot_mat_);
  DEBUG("rot mat " << rot_mat_.str());
  *mobile = rot_mat_.multiply(*mobile);
  DEBUG("mobile after rot " << mobile->str());
  mobile->add(anchor);
  DEBUG("mobile " << mobile->str());

  // update config
  select->get_configuration(system)->update_positions(select->mobile());
}

PerturbConnector::PerturbConnector(std::istream& istr)
  : PerturbMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(3798 == version, "mismatch version: " << version);
}

void PerturbConnector::serialize_perturb_connector_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(3798, ostr);
}

void PerturbConnector::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_connector_(ostr);
}

}  // namespace feasst
