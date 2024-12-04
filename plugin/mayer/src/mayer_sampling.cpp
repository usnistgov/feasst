#include <cmath>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h" // factorial
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "mayer/include/mayer_sampling.h"

namespace feasst {

MayerSampling::MayerSampling(argtype * args) : Criteria(args) {
  class_name_ = "MayerSampling";
  intra_pot_ = integer("intra_potential", args, -1);

  // HWH deprecate
  // Support deprecation warning for old argument name
  if (used("num_trials_per_iteration", *args)) {
    WARN("Metropolis argument num_trials_per_iteration is deprecated. " <<
         "Use trials_per_cycle instead.");
    trials_per_cycle_ = integer("num_trials_per_iteration", args);
  } else {
    trials_per_cycle_ = integer("trials_per_cycle", args, 1e9);
  }
  const int num_beta_taylor = integer("num_beta_taylor", args, 0);
  beta_taylor_.resize(num_beta_taylor);
  for (int ibt = 0; ibt < num_beta_taylor; ++ibt) {
    beta_taylor_[ibt] = *MakeAccumulator({{"num_moments", "2"},
                                          {"max_block_operations", "0"}});
  }
  training_file_ = str("training_file", args, "");
  training_per_write_ = integer("training_per_write", args, 1e4);
}
MayerSampling::MayerSampling(argtype args) : MayerSampling(&args) {
  feasst_check_all_used(args);
}

void MayerSampling::precompute(System * system) {
  system->remove_opt_overlap();
//  const int aniso_index_ = system->configuration().model_params().index("anisotropic");
//  DEBUG("aniso_index " << aniso_index_);
}

bool MayerSampling::is_accepted(
    const System& system,
    Acceptance * acceptance,
    Random * random) {
  check_num_cycles_(trials_per_cycle_);
  ASSERT(system.num_configurations() == 1, "assumes 1 config");
  double energy_new = acceptance->energy_new();
  if (!training_file_.empty()) {
    const Site& mobile = system.configuration().particle(1).site(0);
    Position spherical = mobile.position().spherical();
    DEBUG(mobile.euler().str());
    data_.push_back(std::vector<double>({
      spherical.coord(0), spherical.coord(1), spherical.coord(2),
      mobile.euler().phi(), mobile.euler().theta(), mobile.euler().psi(),
      energy_new}));
    //DEBUG(feasst_str(data_.back()));
    if (static_cast<int>(data_.size()) >= training_per_write_) {
      std::ofstream file;
      file.open(training_file_, std::ofstream::out | std::ofstream::app);
      for (const std::vector<double>& dat : data_) {
        for (const double d : dat) {
          file << d << ",";
          //file << MAX_PRECISION << d << ",";
        }
        file << std::endl;
      }
      file.close();
      data_.clear();
    }
  }
  const double beta = system.thermo_params().beta();
  TRACE("*** MayerSampling ***");
  if (intra_pot_ != -1) {
    TRACE("intra_pot " << intra_pot_);
    const double energy_intra = acceptance->energy_profile_new()[intra_pot_];
    TRACE("energy_intra " << energy_intra);
    if (random->uniform() < std::exp(-beta*energy_intra)) {
      energy_new -= energy_intra;
    } else {
      was_accepted_ = false;
      TRACE("rejected at intra_potential step");
      return was_accepted_;
    }
  }
  const double f12 = std::exp(-beta*energy_new) - 1.;
  TRACE("energy new " << energy_new);
  TRACE("f12 " << f12);
  TRACE("f12old " << f12old_);
  TRACE("acceptance " << std::abs(f12)/std::abs(f12old_));

  if (!acceptance->reject() &&
      (random->uniform() < std::abs(f12)/std::abs(f12old_))) {
    ASSERT(energy_new != 0, "error");
    set_current_energy(acceptance->energy_new());
    set_current_energy_profile(acceptance->energy_profile_new());
    f12old_ = f12;
    was_accepted_ = true;
    TRACE("computing ref");
    f12ref_ = std::exp(-beta*acceptance->energy_ref()) - 1.;
    TRACE("f12ref " << f12ref_);
  } else {
    was_accepted_ = false;
  }
  if (f12old_ < 0) {
    mayer_.accumulate(-1.);
  } else {
    mayer_.accumulate(1.);
  }
  const double mayer_pi = std::abs(f12old_);
  mayer_ref_.accumulate(f12ref_/mayer_pi);
  double energy_old = current_energy();
  if (intra_pot_ != -1) {
    energy_old -= current_energy_profile()[intra_pot_];
  }
  double unebu = std::exp(-beta*energy_old)/mayer_pi;
  for (int ibd = 0; ibd < static_cast<int>(beta_taylor_.size()); ++ibd) {
    unebu *= -energy_old;
    beta_taylor_[ibd].accumulate(unebu/factorial(ibd+1));
  }
  DEBUG("is accepted? " << was_accepted_);
  return was_accepted_;
}

FEASST_MAPPER(MayerSampling,);

void MayerSampling::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(3253, ostr);
  feasst_serialize(f12old_, ostr);
  feasst_serialize(f12ref_, ostr);
  feasst_serialize(trials_per_cycle_, ostr);
  feasst_serialize_fstobj(mayer_, ostr);
  feasst_serialize_fstobj(mayer_ref_, ostr);
  feasst_serialize(intra_pot_, ostr);
  feasst_serialize_fstobj(beta_taylor_, ostr);
  feasst_serialize(training_file_, ostr);
  feasst_serialize(training_per_write_, ostr);
}

MayerSampling::MayerSampling(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3251 && version <= 3253, "unrecognized verison: " << version);
  feasst_deserialize(&f12old_, istr);
  feasst_deserialize(&f12ref_, istr);
  feasst_deserialize(&trials_per_cycle_, istr);
  feasst_deserialize_fstobj(&mayer_, istr);
  feasst_deserialize_fstobj(&mayer_ref_, istr);
  feasst_deserialize(&intra_pot_, istr);
  if (version >= 3252) {
    feasst_deserialize_fstobj(&beta_taylor_, istr);
  }
  if (version >= 3253) {
    feasst_deserialize(&training_file_, istr);
    feasst_deserialize(&training_per_write_, istr);
  }
}

double MayerSampling::second_virial_ratio() const {
  return mayer_.average()/mayer_ref_.average();
}

//MayerSampling::MayerSampling(const Criteria& criteria) {
//  std::stringstream ss;
//  criteria.serialize(ss);
//  *this = MayerSampling(ss);
//}

double MayerSampling::second_virial_ratio_block_stdev() const {
  return std::sqrt(
    pow(mayer_.block_stdev()/mayer_ref_.average(), 2) +
    pow(second_virial_ratio()*mayer_ref_.block_stdev()/mayer_ref_.average(), 2));
}

std::string MayerSampling::write() const {
  std::stringstream ss;
  ss << "{\"second_virial_ratio\": " << second_virial_ratio() << ", "
     << "\"second_virial_ratio_block_stdev\": " << second_virial_ratio_block_stdev();
  if (num_beta_taylor() > 0) {
    ss << ",\"beta_taylor\": [";
    for (int ibt = 0; ibt < num_beta_taylor() + 1; ++ibt) {
      ss << beta_taylor(ibt) << ",";
    }
    ss << "],";
  }
  ss << "}" << std::endl;
  ss << "mayer" << std::endl << mayer().str() << std::endl;
  ss << "mayer_ref" << std::endl << mayer_ref().str() << std::endl;
  return ss.str();
}

double MayerSampling::beta_taylor(const int deriv) const {
  if (deriv == 0) {
    return mayer().average()/mayer_ref().average();
  } else {
    return beta_taylor_[deriv -1].average()/mayer_ref().average();
  }
}

}  // namespace feasst
