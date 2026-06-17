#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "shape/include/sphere.h"
#include "configuration/include/site.h"
#include "configuration/include/configuration.h"
#include "confinement/include/model_table_cartesian_1d.h" // precompute_table1D
#include "confinement/include/model_table_spherical.h"

namespace feasst {

FEASST_MAPPER(ModelTableSphere1D,);

ModelTableSphere1D::ModelTableSphere1D(argtype * args) {
  class_name_ = "ModelTableSphere1D";
  sphere_ = std::make_unique<Sphere>(args);
  table_file_ = str("table_file", args, "");
}
ModelTableSphere1D::ModelTableSphere1D(argtype args) : ModelTableSphere1D(&args) {
  feasst_check_all_used(args);
}
ModelTableSphere1D::~ModelTableSphere1D() {}

void ModelTableSphere1D::precompute(Configuration * config) {
  precompute_table1D(*config, table_file_, &tables_);
}

double ModelTableSphere1D::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  TRACE(tables_.size());
  TRACE(site.type());
  const Table1D* table = tables_[site.type()].get();
  if (!table) {
    return 0.;
  }
  const double radius = sphere_->radius();
  TRACE("radius:" << radius);
  // negative distance means inside shape (cavity)
  const double distance = -sphere_->nearest_distance(wrapped_site);
  TRACE("distance:" << distance);
  if (distance <= NEAR_ZERO) {
    TRACE(MAX_PRECISION << distance);
    TRACE("en: "  << NEAR_INFINITY);
    return NEAR_INFINITY;
  } else if (distance <= radius) {
    const double en = table->forward_difference_interpolation(distance/radius);
    TRACE("en: " << en);
    return en;
  } else if (distance < radius + NEAR_ZERO) {
    const double en = table->forward_difference_interpolation(1.);
    TRACE("en: " << en);
    return en;
  }
  TRACE("en: 0");
  return 0.;
}

void ModelTableSphere1D::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(2675, ostr);
  feasst_serialize(sphere_, ostr);
  feasst_serialize(tables_, ostr);
}

ModelTableSphere1D::ModelTableSphere1D(std::istream& istr) : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2675 && version <= 2675, "unrecognized verison: " << version);
  feasst_deserialize(sphere_, istr);
  feasst_deserialize(&tables_, istr);
}

}  // namespace feasst
