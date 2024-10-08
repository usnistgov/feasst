#include <cmath>  // isnan, pow
#include <string>
#include <fstream>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/utils_math.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "aniso/include/visit_model_inner_nn.h"

namespace feasst {

VisitModelInnerNN::VisitModelInnerNN(argtype * args) : VisitModelInnerTable(args) {
  class_name_ = "VisitModelInnerNN";
}
VisitModelInnerNN::VisitModelInnerNN(argtype args) : VisitModelInnerNN(&args) {
  feasst_check_all_used(args);
}

void VisitModelInnerNN::read_table(const std::string file_name,
    const bool ignore_energy,
    Configuration * config) {
  // HWH implement reading of NN here
}

double VisitModelInnerNN::compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const double e1, const double e2, const double e3, const Configuration& config) const {
  // HWH implement NN predictions here
  double en = 0.;
  return en;
}

FEASST_MAPPER(VisitModelInnerNN,);

VisitModelInnerNN::VisitModelInnerNN(std::istream& istr) : VisitModelInnerTable(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2376, "unrecognized version: " << version);
  // HWH serialize NN here
}

void VisitModelInnerNN::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_table_(ostr);
  feasst_serialize_version(2376, ostr);
  // HWH deserialize NN here
}

}  // namespace feasst
