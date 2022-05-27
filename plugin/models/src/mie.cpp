#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "models/include/mie.h"

namespace feasst {

Mie::Mie(argtype * args) {
  class_name_ = "Mie";
  n_ = dble("n", args, 12);
  m_ = dble("m", args, 6);
  prefactor_ = (n_/(n_ - m_))*std::pow(n_/m_, m_/(n_ - m_));
}
Mie::Mie(argtype args) : Mie(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapMie {
 public:
  MapMie() {
    Mie().deserialize_map()["Mie"] = MakeMie();
  }
};

static MapMie map_mie_ = MapMie();

void Mie::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_mie_(ostr);
}

void Mie::serialize_mie_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(2905, ostr);
  feasst_serialize(n_, ostr);
  feasst_serialize(m_, ostr);
  feasst_serialize(prefactor_, ostr);
}

Mie::Mie(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2905 == version, version);
  feasst_deserialize(&n_, istr);
  feasst_deserialize(&m_, istr);
  feasst_deserialize(&prefactor_, istr);
}

double Mie::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double s_r = sigma/std::sqrt(squared_distance);
  return prefactor_*epsilon*(std::pow(s_r, n_) - std::pow(s_r, m_));
}

}  // namespace feasst
