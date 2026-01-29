#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "shape/include/cylinder.h"
#include "configuration/include/site.h"
#include "configuration/include/configuration.h"
#include "confinement/include/model_table_cartesian.h" // precompute_table1D
#include "confinement/include/model_table_cylindrical.h"

namespace feasst {

FEASST_MAPPER(ModelTableCylinder1D, argtype({{"radius", "1"},
  {"first_point", "0,0,0"}, {"second_point", "0,0,0"}}));

ModelTableCylinder1D::ModelTableCylinder1D(argtype * args) {
  class_name_ = "ModelTableCylinder1D";
  cylinder_ = std::make_unique<Cylinder>(args);
  table_file_ = str("table_file", args, "");
}
ModelTableCylinder1D::ModelTableCylinder1D(argtype args) : ModelTableCylinder1D(&args) {
  feasst_check_all_used(args);
}
ModelTableCylinder1D::~ModelTableCylinder1D() {}

void ModelTableCylinder1D::precompute(Configuration * config) {
  precompute_table1D(*config, table_file_, &tables_);
}

double ModelTableCylinder1D::energy(
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
  const double radius = cylinder_->radius();
  TRACE("radius:" << radius);
  // negative distance means inside shape (cavity)
  const double distance = -cylinder_->nearest_distance(wrapped_site);
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

void ModelTableCylinder1D::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6324, ostr);
  feasst_serialize(cylinder_, ostr);
  feasst_serialize(tables_, ostr);
}

ModelTableCylinder1D::ModelTableCylinder1D(std::istream& istr) : ModelOneBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6324 && version <= 6324, "unrecognized verison: " << version);
  feasst_deserialize(cylinder_, istr);
  feasst_deserialize(&tables_, istr);
}

}  // namespace feasst
