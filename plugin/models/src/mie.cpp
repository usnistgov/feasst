#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "models/include/mie.h"

namespace feasst {

Mie::Mie(const argtype& args) {
  class_name_ = "Mie";
  Arguments args_(args);
  n_ = args_.key("n").dflt("12").dble();
  m_ = args_.key("m").dflt("6").dble();
  prefactor_ = (n_/(n_ - m_))*std::pow(n_/m_, m_/(n_ - m_));
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
  const double sigma = model_params.mixed_sigma()[type1][type2];
  const double epsilon = model_params.mixed_epsilon()[type1][type2];
  const double s_r = sigma/std::sqrt(squared_distance);
  return prefactor_*epsilon*(std::pow(s_r, n_) - std::pow(s_r, m_));
}

}  // namespace feasst
