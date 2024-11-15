#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/extensive_moments.h"

namespace feasst {

FEASST_MAPPER(ExtensiveMoments,);

ExtensiveMoments::ExtensiveMoments(argtype * args) : Analyze(args) {
  max_order_ = integer("max_order", args, 3);
}
ExtensiveMoments::ExtensiveMoments(argtype args) : ExtensiveMoments(&args) {
  feasst_check_all_used(args);
}

void ExtensiveMoments::initialize(MonteCarlo * mc) {
  const int num_ptypes = configuration(mc->system()).num_particle_types();
  resize(max_order_ + 1,
         max_order_ + 1,
         num_ptypes,
         max_order_ + 1,
         num_ptypes,
         &moments_);
  u_p_.resize(max_order_ + 1);
  resize(num_ptypes, max_order_ + 1, &n_i_j_);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string ExtensiveMoments::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  //ss << "max_order," << max_order_ << std::endl;
  return ss.str();
}

void ExtensiveMoments::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  const Criteria& criteria = mc.criteria();
  const int num_ptypes = system.configuration().num_particle_types();

  // precompute powers
  const double energy = criteria.current_energy();
  u_p_[0] = 1.;
  for (int order = 1; order <= max_order_; ++order) {
    u_p_[order] = energy*u_p_[order - 1];
  }
  for (int ptype = 0; ptype < num_ptypes; ++ptype) {
    const int num = system.configuration().num_particles_of_type(ptype);
    n_i_j_[ptype][0] = 1.;
    for (int order = 1; order <= max_order_; ++order) {
      n_i_j_[ptype][order] = num*n_i_j_[ptype][order - 1];
    }
  }

  // accumulate moments
  for (int p = 0; p <= max_order_; ++p) {
  for (int m = 0; m <= max_order_; ++m) {
  for (int k = 0; k < num_ptypes; ++k) {
  for (int j = 0; j <= max_order_; ++j) {
  for (int i = 0; i < num_ptypes; ++i) {
    moments_[p][m][k][j][i].accumulate(n_i_j_[i][j]*n_i_j_[k][m]*u_p_[p]);
  }}}}}
}

std::string ExtensiveMoments::write(const MonteCarlo& mc) {
  std::stringstream ss;
  feasst_serialize_fstobj(moments_, ss);
  return ss.str();
}

void ExtensiveMoments::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1647, ostr);
  feasst_serialize(max_order_, ostr);
  feasst_serialize_fstobj(moments_, ostr);
  feasst_serialize(u_p_, ostr);
  feasst_serialize(n_i_j_, ostr);
}

ExtensiveMoments::ExtensiveMoments(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1647, "mismatch version:" << version);
  feasst_deserialize(&max_order_, istr);
  feasst_deserialize_fstobj(&moments_, istr);
  feasst_deserialize(&u_p_, istr);
  feasst_deserialize(&n_i_j_, istr);
}

ExtensiveMoments::ExtensiveMoments(const Analyze& extensive_moments) {
  std::stringstream ss;
  extensive_moments.serialize(ss);
  *this = ExtensiveMoments(ss);
}

}  // namespace feasst
