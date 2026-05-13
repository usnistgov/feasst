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
  particle_type_names_ = split(str("particle_types", args, ""), ',');
  triplet_ = boolean("triplet", args, false);
}
ExtensiveMoments::ExtensiveMoments(argtype args) : ExtensiveMoments(&args) {
  feasst_check_all_used(args);
}

void ExtensiveMoments::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  const Configuration& config = configuration(mc->system());
  if (particle_type_names_.empty()) {
    // add all particle types
    for (int pt = 0; pt < config.num_particle_types(); ++pt) {
      particle_type_names_.push_back(config.particle_type_to_name(pt));
    }
  }
  for (const std::string& ptn : particle_type_names_) {
    particle_types_.push_back(config.particle_name_to_type(ptn));
  }
  const int num_ptypes = static_cast<int>(particle_types_.size());
  if (triplet_) {
    resize(max_order_ + 1,
           max_order_ + 1,
           num_ptypes,
           max_order_ + 1,
           num_ptypes,
           max_order_ + 1,
           num_ptypes,
           &moments3_);
  } else {
    resize(max_order_ + 1,
           max_order_ + 1,
           num_ptypes,
           max_order_ + 1,
           num_ptypes,
           &moments_);
  }
  u_p_.resize(max_order_ + 1);
  resize(num_ptypes, max_order_ + 1, &n_i_j_);
  printer(header(*mc), output_file(mc->criteria()));
}

std::string ExtensiveMoments::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  if (triplet_) {
    ss << "#{\"max_order\":" << max_order_ << "}" << std::endl
       << "i,j,k,l,m,o,p,average,block_stdev" << std::endl;
  } else {
    ss << "#{\"max_order\":" << max_order_ << "}" << std::endl
       << "i,j,k,m,p,average,block_stdev" << std::endl;
  }
  return ss.str();
}

void ExtensiveMoments::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  const Criteria& criteria = mc.criteria();

  // precompute powers
  const double energy = criteria.current_energy();
  u_p_[0] = 1.;
  for (int order = 1; order <= max_order_; ++order) {
    u_p_[order] = energy*u_p_[order - 1];
  }
  const int num_ptypes = static_cast<int>(particle_types_.size());
  for (int pit = 0; pit < num_ptypes; ++pit) {
    const int ptype = particle_types_[pit];
    const int num = system.configuration().num_particles_of_type(ptype);
    n_i_j_[pit][0] = 1.;
    for (int order = 1; order <= max_order_; ++order) {
      n_i_j_[pit][order] = num*n_i_j_[pit][order - 1];
    }
  }

  // accumulate moments
  if (triplet_) {
    for (int p = 0; p <= max_order_; ++p) {
    for (int o = 0; o <= max_order_; ++o) {
    for (int m = 0; m < num_ptypes; ++m) {
    for (int l = 0; l <= max_order_; ++l) {
    for (int k = 0; k < num_ptypes; ++k) {
    for (int j = 0; j <= max_order_; ++j) {
    for (int i = 0; i < num_ptypes; ++i) {
      moments3_[p][o][m][l][k][j][i].accumulate(n_i_j_[i][j]*n_i_j_[k][l]*n_i_j_[m][o]*u_p_[p]);
    }}}}}}}
  } else {
    for (int p = 0; p <= max_order_; ++p) {
    for (int m = 0; m <= max_order_; ++m) {
    for (int k = 0; k < num_ptypes; ++k) {
    for (int j = 0; j <= max_order_; ++j) {
    for (int i = 0; i < num_ptypes; ++i) {
      moments_[p][m][k][j][i].accumulate(n_i_j_[i][j]*n_i_j_[k][m]*u_p_[p]);
    }}}}}
  }
}

std::string ExtensiveMoments::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  const Configuration& config = configuration(mc.system());
  const int num_ptypes = static_cast<int>(particle_types_.size());
  if (triplet_) {
    for (int p = 0; p <= max_order_; ++p) {
    for (int o = 0; o <= max_order_; ++o) {
    for (int m = 0; m < num_ptypes; ++m) {
    for (int l = 0; l <= max_order_; ++l) {
    for (int k = 0; k < num_ptypes; ++k) {
    for (int j = 0; j <= max_order_; ++j) {
    for (int i = 0; i < num_ptypes; ++i) {
      const Accumulator& acc = moments3_[p][o][m][l][k][j][i];
      ss << i << "," << j << "," << k << "," << l << "," << m << "," << o << ","
         << p << "," << acc.average() << "," << acc.block_stdev() << std::endl;
    }}}}}}}
  } else {
    for (int p = 0; p <= max_order_; ++p) {
    for (int m = 0; m <= max_order_; ++m) {
    for (int k = 0; k < num_ptypes; ++k) {
    for (int j = 0; j <= max_order_; ++j) {
    for (int i = 0; i < num_ptypes; ++i) {
      const Accumulator& acc = moments_[p][m][k][j][i];
      ss << i << "," << j << "," << k << "," << m << "," << p << ","
         << acc.average() << "," << acc.block_stdev() << std::endl;
    }}}}}
  }
  return ss.str();
}

void ExtensiveMoments::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1648, ostr);
  feasst_serialize(max_order_, ostr);
  feasst_serialize_fstobj(moments_, ostr);
  feasst_serialize_fstobj(moments3_, ostr);
  feasst_serialize(u_p_, ostr);
  feasst_serialize(n_i_j_, ostr);
  feasst_serialize(particle_type_names_, ostr);
  feasst_serialize(particle_types_, ostr);
  feasst_serialize(triplet_, ostr);
}

ExtensiveMoments::ExtensiveMoments(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 1647 && version <= 1648, "mismatch version:" << version);
  feasst_deserialize(&max_order_, istr);
  feasst_deserialize_fstobj(&moments_, istr);
  if (version >= 1648) {
    feasst_deserialize_fstobj(&moments_, istr);
  }
  feasst_deserialize(&u_p_, istr);
  feasst_deserialize(&n_i_j_, istr);
  if (version >= 1648) {
    feasst_deserialize(&particle_type_names_, istr);
    feasst_deserialize(&particle_types_, istr);
    feasst_deserialize(&triplet_, istr);
  }
}

ExtensiveMoments::ExtensiveMoments(const Analyze& extensive_moments) {
  std::stringstream ss;
  extensive_moments.serialize(ss);
  *this = ExtensiveMoments(ss);
}

}  // namespace feasst
